### Generating sets of marker SNPs

## This code is designed to find marker SNPs from my prepared inputs of allele counts for archaic and non-archaic haps.
## While it's a bit cumbersome, the approach is very flexible and allows retrieving of a wide range of SNP sets with different properties.


import gzip
import numpy as np
import copy
import bisect
import sys
import os

#sys.path.insert(0, "/home/guy/Dropbox/PythonScripts/")

CWD = os.getcwd()

DAT_DIR = CWD + '/archaic_human_snps/'
AF_DIR = CWD + '/af_freq_files/'
RES_DIR = CWD + '/snp_lists/'

if os.path.isdir(DAT_DIR) == False:
    os.mkdir(DAT_DIR)
if os.path.isdir(AF_DIR) == False:
    os.mkdir(AF_DIR)
if os.path.isdir(RES_DIR) == False:
    os.mkdir(RES_DIR)


def filter_SNPs_on_column_stats(infile_stats, outfile_snps, conditions = []):
    """
    This filters to a subset of SNPs that have certain properties.
    Conditions are [[and],or], ie. [[a,b,c],[d,e]] means either a, b and c are True, or d and e are True.
    The notation of a condition is [[columns], 0/1/2, value] where columns are
    summed and 0-idx, and 0/1/2 indicates below/equal/above. The value can be a
    integer (interpretted as such), a float (interpretted as a proportion,
    of that chunk set), or or a string (in which case it is interpretted
    as a reference to another column).

    Note that I cannot (yet) filter on nearby SNPs etc.
    
    For example, to filter to Deni marker SNPs from a stats file counting
    allelic state in Denisovan chunks, I might want SNPs that are:
    Derived in all Denisovan chunks, ancestral in all human chunks, homo derived
    in the Altai Denisovan and homo ancestral in the Altai Nean.
    [[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[10],1,0],[[11],1,2],[[12],1,2],[[13],1,0]],
    [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0],[[10],1,2],[[11],1,0],[[12],1,0],[[13],1,2]]

    Alternatively, I might prefer not to condition on sharing with the Altai:
    [[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[12],1,2],[[13],1,0]],
    [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0],[[12],1,0],[[13],1,2]]

    Or I could require that an SNP is observed in a chunk >1 rather than >0 times:
    [[[5],1,0],[[6],1,0],[[7],2,1],[[8],2,0],[[9],1,0],[[12],1,2],[[13],1,0]],
    [[[5],1,1],[[6],2,1],[[7],1,0],[[8],1,0],[[9],2,0],[[12],1,0],[[13],1,2]]
    
    """
    open_in = gzip.open if infile_stats[-3:] == '.gz' else open
    open_out = gzip.open if outfile_snps[-3:] == '.gz' else open
    comparative_functions = [lambda a, b: a < b, lambda a, b: a == b, lambda a, b: a > b]
    line_num = 0
    kept_snps = 0
    with open_in(infile_stats, 'rb') as f_in:
        with open_out(outfile_snps, 'wb') as f_out:
            f_out.write('#' + str(conditions) + '\n')
            for line in f_in:
                split_line = line[0:-1].split('\t')
                #print split_line
                if split_line[0] == 'CHR' or split_line[0] == '#CHR' or split_line[0][0] == '#':
                    #header - this isn't a great way to do this, but at least there'll be an error if I try functions on headers...
                    f_out.write(line)
                else:
                    line_num += 1
                    if line_num % 10000 == 1:
                        print 'Read %d SNPs, kept %d at position %s:%s' %(line_num, kept_snps, split_line[0], split_line[1])
                    for condition in conditions:
                        condition_valid = True
                        for requirement in condition:
                            #print requirement
                            if 0 in requirement[0]:
                                objective_value = split_line[requirement[0][0]]
                            else:
                                objective_value = np.sum([float(split_line[i]) for i in requirement[0]])
                            target_value = split_line[int(requirement[2])] if (type(requirement[2]) == str and 0 not in requirement[0]) else requirement[2] if type(requirement[2]) == int else 'proportion'
                            if target_value == 'proportion':
                                #Calculate the proportion.
                                #The column rules are [6,7],[8,9],[10,11],[12,13]
                                if len(requirement[0]) > 0:
                                    raise RuntimeError('Summation of objective_value not allowed if target_value specifies a proportion')
                                objective_column = requirement[0]
                                if objective_column % 2 == 0:
                                    total_observed = int(split_line[objective_column]) + int(split_line[objective_column + 1])
                                else:
                                    total_observed = int(split_line[objective_column]) + int(split_line[objective_column - 1])
                                target_value = float(requirement[2]) * total_observed
                            #print objective_value, target_value, comparative_functions[requirement[1]](objective_value, target_value), comparative_functions[requirement[1]](objective_value, target_value) == False
                            if comparative_functions[requirement[1]](objective_value, target_value) == False:
                                condition_valid = False
                                break
                        #print condition_valid
                        if condition_valid == True:
                            break
                    if condition_valid == True:
                        f_out.write(line)
                        kept_snps += 1
    return None

def filter_SNPs_on_daf(infile_stats, outfile_snps, daf_diag = [], n_resample = 1):
    """
    This filters based on a multiple rejection sampling.
    Construct an ALT allele freq distribution based on altaf_diag.
    Go through each SNP, retaining it if a random number is below the pdf | alt allele freq. (rejection sampling, having re-weighted by the overall alt AF distribution)
    In the archaic conditioned data, I assign a dignostic allele that might be ancestral or derived based on the pattern in human vs denisovan haplotypes, and count the freq of that SNP.
    So, this function builds four lists:
    i) Is the REF or ALT ancestral?
    ii) Is the diagnostic allele REF or ALT?
    These lists have appropriate relative length and Alt AFs according to the sample.
    The important thing is that the Alt AFs 

    For each SNP, I ask which list I attempt to add it to based on the relative number of the four cases. I then ask whether it is actually added based on DAF.
    """
    assert len(outfile_snps) == n_resample
    open_in = gzip.open if infile_stats[-3:] == '.gz' else open
    open_out = gzip.open if outfile_snps[-3:] == '.gz' else open
    n_snps = len(daf_diag)
    #daf_diag.append([alt / float(ref + alt) if ancestral == 0 else ref / float(ref + alt) if ancestral == 1 else np.nan, 0 if int(diag_state) == ancestral else 1, int(diag_state)])
    daf_refancaltdiag = daf_diag[::,0][(daf_diag[::,1] == 0) * (daf_diag[::,2] == 1)] #Alt derived in for an SNP that is 0-ancestral
    daf_altancaltdiag = daf_diag[::,0][(daf_diag[::,1] == 1) * (daf_diag[::,2] == 1)] #Alt ancestral in for an SNP that is 1-ancestral
    daf_refancrefdiag = daf_diag[::,0][(daf_diag[::,1] == 0) * (daf_diag[::,2] == 0)] #Ref ancestral in for an SNP that is 0-ancestral
    daf_altancrefdiag = daf_diag[::,0][(daf_diag[::,1] == 1) * (daf_diag[::,2] == 0)] #Ref derived in for an SNP that is 1-ancestral

    daf_resolution = 50
    daf_refancaltdiag_pdf = np.histogram(daf_refancaltdiag, bins = daf_resolution, range = [0.0,1.0], density = True)
    daf_bins = daf_refancaltdiag_pdf[1][1:]
    daf_refancaltdiag_pdf = daf_refancaltdiag_pdf[0]
    daf_refancaltdiag_pdf = daf_refancaltdiag_pdf / np.sum(daf_refancaltdiag_pdf)
    daf_altancaltdiag_pdf = np.histogram(daf_altancaltdiag, bins = daf_resolution, range = [0.0,1.0], density = True)[0]
    daf_altancaltdiag_pdf = daf_altancaltdiag_pdf / np.sum(daf_altancaltdiag_pdf)
    daf_refancrefdiag_pdf = np.histogram(daf_refancrefdiag, bins = daf_resolution, range = [0.0,1.0], density = True)[0]
    daf_refancrefdiag_pdf = daf_refancrefdiag_pdf / np.sum(daf_refancrefdiag_pdf)
    daf_altancrefdiag_pdf = np.histogram(daf_altancrefdiag, bins = daf_resolution, range = [0.0,1.0], density = True)[0]
    daf_altancrefdiag_pdf = daf_altancrefdiag_pdf / np.sum(daf_altancrefdiag_pdf)

    rel_probs = np.array([len(daf_refancaltdiag), len(daf_altancaltdiag), len(daf_refancrefdiag), len(daf_altancrefdiag)], dtype = float)

    
    rel_probs = rel_probs / float(np.sum(rel_probs))
    print "Relative probability of each category:", rel_probs
    
    rel_probs_cumsum = np.cumsum(rel_probs)

    rel_prob_refanc_0_02 = rel_probs[0] / float(rel_probs[0] + rel_probs[2])
    rel_prob_altanc_1_13 = rel_probs[1] / float(rel_probs[1] + rel_probs[3])

    
    
    daf_list = [daf_refancaltdiag_pdf, daf_altancaltdiag_pdf, daf_refancrefdiag_pdf, daf_altancrefdiag_pdf]

    snp_list = [[], [], [], []]
    

    all_snps = []
    #Read in all data first.
    with open_in(infile_stats, 'rb') as f_in:
        #Decide which list to attempt to add to
        for line in f_in:
            split_line = line[0:-1].split('\t')
            if split_line[0] == 'CHR' or split_line[0] == '#CHR' or split_line[0][0] == '#':
                #header - this isn't a great way to do this, but at least there'll be an error if I try functions on headers...
                pass
            else:
                split_line = line[0:-1].split('\t')
                snp_anc = int(split_line[5])
                if snp_anc != -1:
                    #snp_anc = 0
                    snp_derived = int(split_line[7 - snp_anc]) + int(split_line[9 - snp_anc])
                    snp_ancestral = int(split_line[6 + snp_anc]) + int(split_line[8 + snp_anc])
                    snp_daf = snp_derived / float(snp_derived + snp_ancestral)
                    
                    all_snps.append([split_line[0], split_line[1], snp_anc, snp_daf])

    
    all_snps_dafanc = np.array([[snp[2],snp[3]] for snp in all_snps])
    #I should correct for ref-anc and alt-anc separately as the DAFs are likely to be different
    daf_refanc = np.array([snp[1] for snp in all_snps_dafanc[all_snps_dafanc[::,0] == 0]], dtype = float)
    print "There are %d SNPs in the full list" %(len(all_snps_dafanc))
    #print "The DAFs of SNPs with reference as ancestral are", daf_refanc
    daf_refanc_modifier = np.histogram(daf_refanc, bins = daf_resolution, range = [0.0,1.0], density = True)[0]
    daf_refanc_spectrum = daf_refanc_modifier / np.sum(daf_refanc_modifier)
    print "The observed DAF spectrum for SNPs with reference as ancestral is", daf_refanc_spectrum
    daf_refanc_modifier = 1.0 / (daf_refanc_spectrum)

    daf_altanc = np.array([snp[1] for snp in all_snps_dafanc[all_snps_dafanc[::,0] == 1]], dtype = float)
    #print "The DAFs of SNPs with alt as ancestral are", daf_altanc
    daf_altanc_modifier = np.histogram(daf_altanc, bins = daf_resolution, range = [0.0,1.0], density = True)[0]
    daf_altanc_spectrum = daf_altanc_modifier / np.sum(daf_altanc_modifier)
    print "The observed DAF spectrum for SNPs with alt as ancestral is", daf_altanc_spectrum
    daf_altanc_modifier = 1.0 / (daf_altanc_spectrum)

    daf_list_modified = [daf_refancaltdiag_pdf * daf_refanc_modifier, daf_altancaltdiag_pdf * daf_altanc_modifier, daf_refancrefdiag_pdf * daf_refanc_modifier, daf_altancrefdiag_pdf * daf_altanc_modifier]
    daf_list_modified = [daf_modified / np.max(daf_modified[np.isinf(daf_modified) == False]) for daf_modified in daf_list_modified]

    """
    print "Modified main sampling probs:", daf_list_modified[0]
    
    print "Target refancaltdiag allele freq spectrum:", daf_refancaltdiag_pdf
    print "Expected refancaltdiag allele freq spectrum:", (daf_list_modified[0] * daf_refanc_spectrum) / np.sum((daf_list_modified[0] * daf_refanc_spectrum))

    print "Target altancaltdiag allele freq spectrum:", daf_altancaltdiag_pdf
    print "Expected altancaltdiag allele freq spectrum:", (daf_list_modified[1] * daf_altanc_spectrum) / np.sum((daf_list_modified[1] * daf_altanc_spectrum))

    print "Target refancrefdiag allele freq spectrum:", daf_refancrefdiag_pdf
    print "Expected refancrefdiag allele freq spectrum:", (daf_list_modified[2] * daf_refanc_spectrum) / np.sum((daf_list_modified[2] * daf_refanc_spectrum))

    print "Target altancrefdiag allele freq spectrum:", daf_altancrefdiag_pdf
    print "Expected altancrefdiag allele freq spectrum:", (daf_list_modified[3] * daf_altanc_spectrum) / np.sum((daf_list_modified[3] * daf_altanc_spectrum))
    """
    for snp in all_snps:
        snp_daf_idx = bisect.bisect_left(daf_bins, snp[3])
        #Split into 0-anc and 1-anc.
        if snp[2] == 0:
            #There should be plenty of SNPs to choose from if I just split 50/50, unless perhaps the ascertainment population is very small
            #daf_idx = 0 if np.random.rand() < rel_prob_refanc_0_02 else 2
            daf_idx = 0 if np.random.rand() < 0.5 else 2
        elif snp[2] == 1:
            #There should be plenty of SNPs to choose from if I just split 50/50, unless perhaps the ascertainment population is very small
            #daf_idx = 1 if np.random.rand() < rel_prob_altanc_1_13 else 3
            daf_idx = 1 if np.random.rand() < 0.5 else 3
        #print daf_idx, snp_daf, daf_list[daf_idx][snp_daf_idx]
        if np.random.rand() < daf_list_modified[daf_idx][snp_daf_idx]:
            snp_data = [snp[0], snp[1], '1' if daf_idx in [0,1] else '0', snp[3], snp[2]] ##Chr, pos, diag, DAF, anc
            #print snp_data
            snp_list[daf_idx].append(snp_data)
    
    # Now, choose from the snp_lists! According to the expectations, ideally without replacement, and write out these diagnostic SNPs.
    for resamp in range(n_resample):
        with open_out(outfile_snps[resamp] + '.snps', 'wb') as f_out_snps:
            with open_out(outfile_snps[resamp] + '.snpsDiag', 'wb') as f_out_diag:
                successful_sample = False
                attempts = 0
                while successful_sample == False:
                    #print rel_probs, np.sum(rel_probs)
                    num_samp = np.random.multinomial(len(daf_diag), rel_probs)
                    #print num_samp
                    try:
                        resamp_snps = []
                        for i in range(len(rel_probs)):
                            if attempts < 10:
                                choices = np.random.choice(np.arange(len(snp_list[i])), num_samp[i], replace = False)
                            else:
                                print "Failed after 10 attempts, sampling with replacement"
                                choices = np.random.choice(np.arange(len(snp_list[i])), num_samp[i], replace = True)
                            resamp_snps.extend([snp_list[i][j] for j in choices])
                            #resamp_snps.extend(list(np.random.choice(snp_list[i], num_samp[i])))
                        resamp_chr_pos_snps = np.array([[int(chr_pos_snp[0].replace('chr', '')), int(chr_pos_snp[1])] for chr_pos_snp in resamp_snps])
                        resamp_order = np.lexsort((resamp_chr_pos_snps[::,1], resamp_chr_pos_snps[::,0]))
                        resamp_dafs = []
                        resamp_diagAFs = []
                        
                        for resamp_snp_idx in resamp_order:
                            f_out_snps.write('\t'.join(resamp_snps[resamp_snp_idx][0:2]) + '\n') #Chr pos
                            f_out_diag.write('\t'.join(resamp_snps[resamp_snp_idx][0:3]) + '\n') #Chr pos diag
                            resamp_dafs.append(resamp_snps[resamp_snp_idx][3])
                            resamp_diagAFs.append(resamp_snps[resamp_snp_idx][3] if int(resamp_snps[resamp_snp_idx][2]) + resamp_snps[resamp_snp_idx][4] == 1 else 1.0 - resamp_snps[resamp_snp_idx][3])
                            
                        print "Average DAF of resampled SNPs in target population is %.3f compared to %.3f expected." %(np.average(resamp_dafs), np.average(daf_diag[::,0]))
                        print "Average Diagnostic Freq of resampled SNPs in target population is %.3f compared to %.3f expected." %(np.average(resamp_diagAFs), np.average(np.concatenate((daf_diag[::,0][(daf_diag[::,1] + daf_diag[::,2]) == 1], 1.0 - daf_diag[::,0][(daf_diag[::,1] + daf_diag[::,2]) != 1]))))
                        successful_sample = True
                    except ValueError:
                        # I've seen this happen when the true Deni alleles are ancestral REF alleles re-entering the population where the derived ALT allele is at high frequency.
                        # For now let's take this as a failure and sample with replacement after a few goes.
                        print "Failed at %d-th attempt" %(attempts)
                        print num_samp
                        print [len(snp_list[i]) for i in range(4)]
                        print rel_probs
                        attempts += 1
                #pl.hist(daf_refancaltdiag, bins = 50, range = [0.0,1.0], density = True, label = 'expected')
                #pl.show()
                #pl.hist(tmp[0], bins = 50, range = [0.0,1.0], density = True, label = 'resample')
                #pl.show()

    return resamp_snps

def winnow_snp_list(infile_bed, outfile, min_dist = 10000, conditions = None, exclusion_dict = {}):
    """
    This requires a minimum distance between SNPs, and optionally for a chr_pos SNP string not to be in a dictionary of banned SNPs (usually, to remove SNPs that are polymorphic in Africa)
    
    The default method is to iteratively check the distance between the current
    SNP[i] and the SNP[i+1] and and tag SNP[i+1] for removal if the distnace < x
    before checkign [i+2,...+k] etc. until an SNP isn't deleted. Then move to
    position [i+k] and continue.
    
    Alternatively, I can preferentially choose an SNP based on some operation,
    For example, the highest overall Deni introgressed frequency.
    Draws are solved by the default method.

    Example condition for maximum total Deni:
    [1,[6,7]] #The 1 indicated maximum and the [6,7] is the summation of columns

    It may not matter all that much how SNPs are removed, but I do want the
    output to be consistent.
    """
    open_in = gzip.open if infile_bed[-3:] == '.gz' else open
    open_out = gzip.open if outfile[-3:] == '.gz' else open
    if min_dist == 0:
        with open_in(infile_bed, 'rb') as f_in:
            with open_out(outfile, 'wb') as f_out:
                for line in f_in:
                    if line[0] == '#':
                        f_out.write(line)
                    else:
                        split_line = line[0:-1].split('\t')
                        chr_pos = '_'.join([split_line[0].replace('chr', ''), split_line[1]])
                        if chr_pos not in exclusion_dict:
                            f_out.write(line)
    else:
        headers = []
        winnowed = []
        line_choice = [None, None]
        tot_snps = 0
        tot_dropped = 0
        with open_in(infile_bed, 'rb') as f_in:
            curr_chrom = None
            curr_pos = None
            next_pos = None
            next_chrom = None
            for line in f_in:
                split_line = line[0:-1].split('\t')
                if line[0] == '#' or split_line[0] == 'CHR':
                    headers.append(line)
                else:
                    tot_snps += 1
                    #print curr_pos, next_pos
                    #print split_line
                    #print curr_chrom is None
                    chr_pos = '_'.join([split_line[0].replace('chr', ''), split_line[1]])
                    if chr_pos in exclusion_dict:
                        #Skip this line, the SNP is in the exclusion dictionary
                        pass
                    else:
                        if (curr_chrom is None) or split_line[0] != curr_chrom:
                            #New chromosome!
                            #Add any remainins SNPs
                            if line_choice[0] is not None:
                                winnowed.append(line_choice[0])
                            if line_choice[1] is not None:
                                winnowed.append(line_choice[1])
                            #Then reset
                            curr_chrom = split_line[0]
                            curr_pos = int(split_line[1])
                            next_pos = None
                            next_chrom = None
                            line_choice = [split_line, None]
                        elif line_choice[1] is None:
                            #Update the choice
                            next_chrom = split_line[0]
                            next_pos = int(split_line[1])
                            line_choice[1] = split_line
                        #print line_choice
                        #print curr_pos, next_pos
                        if (next_pos is not None) and (next_pos - curr_pos) > min_dist:
                            #print next_pos - curr_pos
                            #SNPs are far enough apart
                            winnowed.append(line_choice[0])
                            line_choice[0] = copy.deepcopy(line_choice[1])
                            line_choice[1] = None
                            curr_chrom = copy.deepcopy(next_chrom)
                            curr_pos = copy.deepcopy(next_pos)
                            next_pos = None
                            next_chrom = None
                            #print "Added"
                            #print winnowed[-1]
                            #print line_choice
                        elif next_pos is not None:
                            #Make my choice!
                            #print next_pos - curr_pos
                            #print "Too close"
                            tot_dropped += 1
                            if type(conditions) is type(None):
                                line_choice[1] = None
                                next_pos = None
                                next_chrom = None
                            else:
                                curr_val = np.sum([float(line_choice[0][i]) for i in conditions[1]])
                                next_val = np.sum([float(line_choice[1][i]) for i in conditions[1]])
                                #print curr_val, next_val
                                if curr_val == next_val:
                                    line_choice[1] = None
                                    next_pos = None
                                    next_chrom = None
                                elif (curr_val > next_val and conditions[0] == 1) or (curr_val < next_val and conditions[0] == 0):
                                    #keeping highest which is current or lowest which is current
                                    line_choice[1] = None
                                    next_pos = None
                                    next_chrom = None
                                elif (curr_val > next_val and conditions[0] == 0) or (curr_val < next_val and conditions[0] == 1):
                                    #keeping highest which is next or lowest which is next
                                    line_choice[0] = copy.deepcopy(line_choice[1])
                                    curr_pos = copy.deepcopy(next_pos)
                                    curr_chrom = copy.deepcopy(next_chrom)
                                    line_choice[1] = None
                                    next_pos = None
                                    next_chrom = None
                            #print line_choice
                        else:
                            pass
        if line_choice[0] is not None:
            winnowed.append(line_choice[0])
        with open_out(outfile, 'wb') as f_out:
            for line in headers:
                f_out.write(line)
            for line in winnowed:
                f_out.write('\t'.join(line) + '\n')
        print "Dropped %d of %d SNPs with minimum distance of %d" %(tot_dropped, tot_snps, min_dist)
    return None
    
def bed_to_snp_list(infile_bed, outfile = None, diagnostic_col_over_col = [7,6]):
    """
    This creates a list of SNPs (chr\tpos) to use as the --regions_file in bcftools.
    It also creates a list of SNPs (chr\tpos\taltref) for interpretting the frequency output.
    For example, if the ALT is the diagnostic allele for an SNP it is important to note this.
    
    altref is determined using diagnostic_col_over_col. The argument is that I am choosing SNPs
    based on one column being greater than the other. For example, if 0-idx column 7 is the
    count of the alt state in the tagged haplotypes and column 6 is the count of the ref state
    then diagnostic_col_over_col = [7,6] says that the alt is diagnostic if 7 > 6 (usual case)
    and the ref is diagnostic otherwise.
    """
    open_in = gzip.open if infile_bed[-3:] == '.gz' else open
    if type(outfile) != type(None):
        open_out = gzip.open if outfile[-3:] == '.gz' else open
    else:
        open_out = None
    to_write = []
    daf_diag = []
    with open_in(infile_bed, 'rb') as f_in:
        for line in f_in:
            split_line = line[0:-1].split('\t')
            if line[0] == '#' or split_line[0] == 'CHR':
                pass
            else:
                diag_state = '1' if int(split_line[diagnostic_col_over_col[0]]) > int(split_line[diagnostic_col_over_col[1]]) else '0'
                to_write.append(list(split_line[0:2]) + [diag_state])
                ancestral = int(split_line[5])
                ref = int(split_line[6]) + int(split_line[8])
                alt = int(split_line[7]) + int(split_line[9])
                if ancestral != -1:
                    #This is Alt allele freq
                    #altaf_diag.append([alt / float(ref + alt), ancestral, int(diag_state)])
                    daf_diag.append([alt / float(ref + alt) if ancestral == 0 else ref / float(ref + alt) if ancestral == 1 else np.nan, ancestral, int(diag_state)])
                    #daf_diag.append([alt / float(ref + alt) if ancestral == 0 else ref / float(ref + alt) if ancestral == 1 else np.nan, 0 if int(diag_state) == ancestral else 1, int(diag_state)])
    if type(outfile) != type(None):
        with open_out(outfile[:-3] + '.snps.gz' if outfile[:-3] == '.gz' else outfile + '.snps', 'wb') as f_out:
            for line in to_write:
                f_out.write('\t'.join(line[0:2]) + '\n')
        with open_out(outfile[:-3] + '.snpsDiag.gz' if outfile[:-3] == '.gz' else outfile + '.snpsDiag', 'wb') as f_out:
            for line in to_write:
                #These are the actual SNPs used, after all winnowing, conditioning etc. So I should base my resample analyses on SNPs with similar DAF to these ones.
                f_out.write('\t'.join(line) + '\n')
    return np.array(daf_diag)


def generate_deni_diagnostic(ascertainment_population = 'subset_papuaExclFrancoisUVBaining', overwrite = False, africa_freq_filebase = AF_DIR + "chr%d.AD.AN_highQ_only_wo_N_biSNP.ANC.479.haps.frequency.subset_africa_SSonly.frq.gz", africa_condition = False, resample_daf = False):
    #I'm just going to profile the capacity of Deni diangostic SNPs (with no African conditioning) to mimic, approximately, the Deni proportions in individuals
    #Firstly, I will condition on having to be the same as the Deni.
    #Secondly, I will not condition on this
    #Thirdly, I will allow some portion of Deni to be missed, i.e. 1 or 2 observations of a variant
    #(I need to implement 4 -> condition on Africa low freq)
    
    test_conditions = [
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[10],1,0],[[11],1,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0],[[10],1,2],[[11],1,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,1],[[8],2,0],[[9],1,0],[[10],1,0],[[11],1,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,1],[[7],1,0],[[8],1,0],[[9],2,0],[[10],1,2],[[11],1,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,1],[[8],2,0],[[9],1,0],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,1],[[7],1,0],[[8],1,0],[[9],2,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],0,2],[[10],1,0],[[11],1,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],0,2],[[9],2,0],[[10],1,2],[[11],1,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],0,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],0,2],[[9],2,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0]]],
        [[[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0]],
         [[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0]]],
        [[[[5],1,0],[[6],0,2],[[7],2,1],[[8],2,0],[[9],0,2],[[10],1,0],[[11],1,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,1],[[7],0,2],[[8],0,2],[[9],2,0],[[10],1,2],[[11],1,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[6],0,2],[[7],2,1],[[8],2,0],[[9],0,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,1],[[7],0,2],[[8],0,2],[[9],2,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[10],1,0],[[11],1,2],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,0],[[8],1,0],[[9],2,0],[[10],1,2],[[11],1,0],[[12],1,0],[[13],1,2]]],
        [[[[5],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[12],1,2],[[13],1,0]],
         [[[5],1,1],[[6],2,0],[[8],1,0],[[9],2,0],[[12],1,0],[[13],1,2]]]
        ]
    test_names = ['IDenDerFix_HumDerMax0_ADenDer_ANeaAnc', # introgressed chunks: fixed derived, human chunks: fixed ancestral, Altai Deni: homozygote derived, Altai Nean: homozygoste ancestral
                  'IDenDerFix_HumDerMax0_ADenAny_ANeaAnc', # introgressed chunks: fixed derived, human chunks: fixed ancestral, Altai Deni: any, Altai Nean: homozygoste ancestral
                  'IDenDerFix2_HumDerMax0_ADenDer_ANeaAnc', # introgressed chunks: fixed derived (> 1 observed), human chunks: fixed ancestral, Altai Deni: homozygote derived, Altai Nean: homozygoste ancestral
                  'IDenDerFix2_HumDerMax0_ADenAny_ANeaAnc', # introgressed chunks: fixed derived (> 1 observed), human chunks: fixed ancestral, Altai Deni: any, Altai Nean: homozygoste ancestral
                  'IDenDerFix_HumDerMax1_ADenDer_ANeaAnc', # introgressed chunks: fixed derived, human chunks: max 1 derived, Altai Deni: homozygote derived, Altai Nean: homozygoste ancestral
                  'IDenDerFix_HumDerMax1_ADenAny_ANeaAnc', # introgressed chunks: fixed derived, human chunks: max 1 derived, Altai Deni: any, Altai Nean: homozygoste ancestral
                  'IDenDerFix_HumDerMax0_ADenAny_ANeaAny', # introgressed chunks: fixed derived, human chunks: fixed ancestral, Altai Deni: any, Altai Nean: any
                  'IDenAlleleFix_HumAlleleMax0_ADenAny_ANeaAny', # introgressed chunks: fixed anc or der, human chunks: fixed alternative to introgressed, Altai Deni: any, Altai Nean: any
                  'IDenDerMin2AncMax1_HumDerMax1_ADenDer_ANeaAnc', # introgressed chunks: min 2 der max 1 anc, human chunks: max 1 derived, Altai Deni: homozygote derived, Altai Nean: homozygote ancestral
                  'IDenDerMin2AncMax1_HumDerMax1_ADenAny_ANeaAnc', # introgressed chunks: min 2 der max 1 anc, human chunks: max 1 derived, Altai Deni: any, Altai Nean: homozygote ancestral
                  'IDenDerMin1_HumDerMax0_ADenDer_ANeaAnc', # introgressed chunks: min 1 der any number anc, human chunks: fixed ancestral, Altai Deni: homozygote derived, Altai Nean: homozygote ancestral
                  'IDenDerMin1_HumDerMax0_ADenAny_ANeaAnc'] # introgressed chunks: min 1 der any number anc, human chunks: fixed ancestral, Altai Deni: any, Altai Nean: homozygote ancestral
    
    assert len(test_names) == len(test_conditions)
    if africa_condition == True:
        #This is used to remove SNPs that are polymorphic in Africa (note; if the SNP is fixed in Africa in the Papuan Deni state then I'd not lose that SNP.)
        africa_snp_freq_dict = {}
        open_africa = gzip.open if africa_freq_filebase[-3:] == '.gz' else open
        for chrom in range(1,23):
            with open_africa(africa_freq_filebase %(chrom), 'rb') as out_in:
                headers = out_in.readline()
                for line in out_in:
                    split_line = line[0:-1].split('\t')
                    snp_id = '_'.join([split_line[0].replace('chr', ''), split_line[1]])
                    snp_freq_ref = float(split_line[4].split(':')[1])
                    snp_freq_alt = float(split_line[5].split(':')[1])
                    if snp_freq_ref != 0.0 and snp_freq_ref != 1.0:
                        africa_snp_freq_dict[snp_id] = None
        print "Read in %d africa polymorphic SNP frequencies" %(len(africa_snp_freq_dict))
    else:
        africa_snp_freq_dict = {}
    
    for winnow in ['0', '5000']:
        for test_idx in range(len(test_names)):
            test_condition = test_conditions[test_idx]
            test_name = test_names[test_idx]
            print "Testing: ", test_name
            #The first step is to filter the SNPs down to my preferred conditions
            print "Filtering SNPs using method %s" %(test_name)
            if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.BED.gz' %(ascertainment_population, test_name)):
                print "Skipping step, filtered SNPs already exists"
            else:
                filter_SNPs_on_column_stats(infile_stats = DAT_DIR + '/SNPstats.%s.deniHCSS35unique.BED.gz' %(ascertainment_population),
                                        outfile_snps = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.BED.gz' %(ascertainment_population, test_name),
                                        conditions = test_condition)

            #There is an optional winnowing step. This can drop SNPs that are too close to reduce LD effects, and also drop SNPs that are in a supplied exclusion dict
            if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.%sw%s.BED.gz' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow)) == True:
                print "Skipping step, winnowing; file already exists"
            else:
                winnow_snp_list(infile_bed = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.BED.gz' %(ascertainment_population, test_name),
                                outfile = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.%sw%s.BED.gz' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow),
                                min_dist = int(winnow),
                                conditions = [1,[6,7]],
                                exclusion_dict = africa_snp_freq_dict)
            
            #The third step is to generate an snp file of locations and the diagnostic state
            print "Generating diagnostic SNP list"
            if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.%sw%s.snps' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow)) and os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.%sw%s.snpsDiag' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow)) and resample_daf == False:
                print "Skipping step, diagnostic SNP list already exists"
            else:
                daf_diag = bed_to_snp_list(infile_bed = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.%sw%s.BED.gz' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow),
                            outfile = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%s.%sw%s' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow) if resample_daf == False else None,
                            diagnostic_col_over_col = [7,6])

            #If I am sampling based on DAF of SNPs and the alt/ref diagnostic allele then generate resample_daf downsampled replicate lists of SNPs.
            outfile_snps_list = []
            resample_daf_required = copy.copy(resample_daf)
            if resample_daf != False:
                print "Constructing resampled SNP sets based on ascertained DAF distribution"
                #1. Winnowing step FIRST this time. This at least means no African polymorphic (if required) and all SNPs > some minimum distance apart.
                if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.BED.gz' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow)):
                    print "Skipping winnowing entire SNP set step, file already exists"
                    pass
                else:
                    print "Winnowing entire SNPs set to suitable conditions (i.e. window size, excluding African polymorphic if required)"
                    winnow_snp_list(infile_bed = DAT_DIR + 'SNPstats.%s.deniHCSS35unique.BED.gz' %(ascertainment_population),
                                outfile = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.BED.gz' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow),
                                min_dist = int(winnow),
                                conditions = None,
                                exclusion_dict = africa_snp_freq_dict)
                #2. Do rejection sampling based on DAF on the winnowed set. Randomly choose subsets of this list and write these with diagnostic alleles indicated appropriately.
                if overwrite is False and np.sum([os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.resampleto-%s.%d.snpsDiag' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, i)) for i in range(resample_daf)]) == resample_daf:
                    print "Skipping SNP resampling as %d resampled sets already exist" %(resample_daf)
                    outfile_snps_list = [RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.resampleto-%s.%d' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_idx) for resample_daf_idx in range(resample_daf)]
                    resample_daf_required = copy.copy(resample_daf)
                else:
                    print "Resampling %d sets of SNPs from full winnowed SNP set" %(resample_daf)
                    if overwrite is False:
                        resample_daf_required  = 0
                        while resample_daf_required < resample_daf:
                            if os.path.exists(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.resampleto-%s.%d.snps' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_required)):
                                pass
                            else:
                                outfile_snps_list.append(RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.resampleto-%s.%d' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_required))
                            resample_daf_required += 1
                        resample_daf_required = len(outfile_snps_list)
                    else:
                        outfile_snps_list = [RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.resampleto-%s.%d' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_idx) for resample_daf_idx in range(resample_daf)]
                        resample_daf_required = copy.copy(resample_daf)
                    print "Constructing %d DAF-matching, Anc/Ref ancestral-matching SNP sets" %(resample_daf_required)
                    daf_diag = filter_SNPs_on_daf(infile_stats = RES_DIR + '/SNPstats.%s.deniHCSS35unique.%sw%s.BED.gz' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow),
                                            outfile_snps = outfile_snps_list,
                                            daf_diag = daf_diag,
                                            n_resample = resample_daf_required)
    return None


def generate_nean_diagnostic(ascertainment_population = 'subset_papuaExclFrancois', overwrite = False, africa_freq_filebase = "/media/guy/sf_Genetic_Data/Indonesia_Diversity175/multicall3_mask1/phase/VCFsWithArchaic/frequency/chr%d.AD.AN_highQ_only_wo_N_biSNP.ANC.479.haps.frequency.subset_africa_SSonly.frq.gz", africa_condition = False, resample_daf = False):
    #I'm just going to profile the capacity of Nean diangostic SNPs (with no African conditioning) to mimic, approximately, the Nean proportions in individuals
    #I'm trying 12 different ascertainment schemes.
    test_conditions = [
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[12],1,0],[[13],1,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0],[[12],1,2],[[13],1,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,1],[[8],2,0],[[9],1,0],[[12],1,0],[[13],1,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,1],[[7],1,0],[[8],1,0],[[9],2,0],[[12],1,2],[[13],1,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,1],[[8],2,0],[[9],1,0],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,1],[[7],1,0],[[8],1,0],[[9],2,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],0,2],[[12],1,0],[[13],1,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],0,2],[[9],2,0],[[12],1,2],[[13],1,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],0,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],0,2],[[9],2,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0]],
         [[[5],1,1],[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0]]],
        [[[[6],1,0],[[7],2,0],[[8],2,0],[[9],1,0]],
         [[[6],2,0],[[7],1,0],[[8],1,0],[[9],2,0]]],
        [[[[5],1,0],[[6],0,2],[[7],2,1],[[8],2,0],[[9],0,2],[[12],1,0],[[13],1,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,1],[[7],0,2],[[8],0,2],[[9],2,0],[[12],1,2],[[13],1,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[6],0,2],[[7],2,1],[[8],2,0],[[9],0,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,1],[[7],0,2],[[8],0,2],[[9],2,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[12],1,0],[[13],1,2],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,0],[[8],1,0],[[9],2,0],[[12],1,2],[[13],1,0],[[10],1,0],[[11],1,2]]],
        [[[[5],1,0],[[7],2,0],[[8],2,0],[[9],1,0],[[10],1,2],[[11],1,0]],
         [[[5],1,1],[[6],2,0],[[8],1,0],[[9],2,0],[[10],1,0],[[11],1,2]]]
        ]
    test_names = ['INeaDerFix_HumDerMax0_ANeaDer_ADenAnc',
                  'INeaDerFix_HumDerMax0_ANeaAny_ADenAnc',
                  'INeaDerFix2_HumDerMax0_ANeaDer_ADenAnc',
                  'INeaDerFix2_HumDerMax0_ANeaAny_ADenAnc',
                  'INeaDerFix_HumDerMax1_ANeaDer_ADenAnc',
                  'INeaDerFix_HumDerMax1_ANeaAny_ADenAnc',
                  'INeaDerFix_HumDerMax0_ANeaAny_ADenAny',
                  'INeaAlleleFix_HumAlleleMax0_ANeaAny_ADenAny',
                  'INeaDerMin2AncMax1_HumDerMax1_ANeaDer_ADenAnc',
                  'INeaDerMin2AncMax1_HumDerMax1_ANeaAny_ADenAnc',
                  'INeaDerMin1_HumDerMax0_ANeaDer_ADenAnc',
                  'INeaDerMin1_HumDerMax0_ANeaAny_ADenAnc']
    assert len(test_names) == len(test_conditions)
    if africa_condition == True:
        #This is used to remove SNPs that are polymorphic in Africa (note; if the SNP is fixed in Africa in the Papuan Deni state then I'd not lose that SNP.)
        africa_snp_freq_dict = {}
        open_africa = gzip.open if africa_freq_filebase[-3:] == '.gz' else open
        for chrom in range(1,23):
            with open_africa(africa_freq_filebase %(chrom), 'rb') as out_in:
                headers = out_in.readline()
                for line in out_in:
                    split_line = line[0:-1].split('\t')
                    snp_id = '_'.join([split_line[0].replace('chr', ''), split_line[1]])
                    snp_freq_ref = float(split_line[4].split(':')[1])
                    snp_freq_alt = float(split_line[5].split(':')[1])
                    if snp_freq_ref != 0.0 and snp_freq_ref != 1.0:
                        africa_snp_freq_dict[snp_id] = None
        print "Read in %d africa polymorphic SNP frequencies" %(len(africa_snp_freq_dict))
    else:
        africa_snp_freq_dict = {}

    
    for winnow in ['0', '5000']:
        for test_idx in range(len(test_names)):
            test_condition = test_conditions[test_idx]
            test_name = test_names[test_idx]
            print "Testing: ", test_name
            #The first step is to filter the SNPs down to my preferred conditions
            print "Filtering SNPs using method %s" %(test_name)
            if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.BED.gz' %(ascertainment_population, test_name)):
                print "Skipping step, filtered SNPs already exists"
            else:
                filter_SNPs_on_column_stats(infile_stats = DAT_DIR + '/SNPstats.%s.neanHCSS35unique.BED.gz' %(ascertainment_population),
                                        outfile_snps = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.BED.gz' %(ascertainment_population, test_name),
                                        conditions = test_condition)

            #There is an optional winnowing step. This can drop SNPs that are too close to reduce LD effects, and also drop SNPs that are in a supplied exclusion dict
            if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.%sw%s.BED.gz' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow)) == True:
                print "Skipping step, winnowing; file already exists"
            else:
                winnow_snp_list(infile_bed = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.BED.gz' %(ascertainment_population, test_name),
                                outfile = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.%sw%s.BED.gz' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow),
                                min_dist = int(winnow),
                                conditions = [1,[6,7]],
                                exclusion_dict = africa_snp_freq_dict)
            
            #The third step is to generate an snp file of locations and the diagnostic state
            print "Generating diagnostic SNP list"
            if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.%sw%s.snps' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow)) and os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.%sw%s.snpsDiag' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow)) and resample_daf == False:
                print "Skipping step, diagnostic SNP list already exists"
            else:
                daf_diag = bed_to_snp_list(infile_bed = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.%sw%s.BED.gz' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow),
                            outfile = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%s.%sw%s' %(ascertainment_population, test_name, 'noAfPoly.' if africa_condition == True else '', winnow) if resample_daf == False else None,
                            diagnostic_col_over_col = [7,6])

            #If I am sampling based on DAF of SNPs and the alt/ref diagnostic allele then generate resample_daf downsampled replicate lists of SNPs.
            outfile_snps_list = []
            resample_daf_required = copy.copy(resample_daf)
            if resample_daf != False:
                print "Constructing resampled SNP sets based on ascertained DAF distribution"
                #1. Winnowing step FIRST this time. This at least means no African polymorphic (if required) and all SNPs > some minimum distance apart.
                if overwrite is False and os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.BED.gz' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow)):
                    print "Skipping winnowing entire SNP set step, file already exists"
                    pass
                else:
                    print "Winnowing entire SNPs set to suitable conditions (i.e. window size, excluding African polymorphic if required)"
                    winnow_snp_list(infile_bed = DAT_DIR + '/SNPstats.%s.neanHCSS35unique.BED.gz' %(ascertainment_population),
                                outfile = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.BED.gz' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow),
                                min_dist = int(winnow),
                                conditions = None,
                                exclusion_dict = africa_snp_freq_dict)
                #2. Do rejection sampling based on DAF on the winnowed set. Randomly choose subsets of this list and write these with diagnostic alleles indicated appropriately.
                if overwrite is False and np.sum([os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.resampleto-%s.%d.snpsDiag' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, i)) for i in range(resample_daf)]) == resample_daf:
                    print "Skipping SNP resampling as %d resampled sets already exist" %(resample_daf)
                    outfile_snps_list = [RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.resampleto-%s.%d' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_idx) for resample_daf_idx in range(resample_daf)]
                    resample_daf_required = copy.copy(resample_daf)
                else:
                    print "Resampling %d sets of SNPs from full winnowed SNP set" %(resample_daf)
                    if overwrite is False:
                        resample_daf_required  = 0
                        while resample_daf_required < resample_daf:
                            if os.path.exists(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.resampleto-%s.%d.snps' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_required)):
                                pass
                            else:
                                outfile_snps_list.append(RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.resampleto-%s.%d' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_required))
                            resample_daf_required += 1
                        resample_daf_required = len(outfile_snps_list)
                    else:
                        outfile_snps_list = [RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.resampleto-%s.%d' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow, test_name, resample_daf_idx) for resample_daf_idx in range(resample_daf)]
                        resample_daf_required = copy.copy(resample_daf)
                    print "Constructing %d DAF-matching, Anc/Ref ancestral-matching SNP sets" %(resample_daf_required)
                    daf_diag = filter_SNPs_on_daf(infile_stats = RES_DIR + '/SNPstats.%s.neanHCSS35unique.%sw%s.BED.gz' %(ascertainment_population, 'noAfPoly.' if africa_condition == True else '', winnow),
                                            outfile_snps = outfile_snps_list,
                                            daf_diag = daf_diag,
                                            n_resample = resample_daf_required)
    return None


if True:
    for ascertainment_population in ['continent_sasia', 'continent_easia', 'subset_papuaExclFrancoisUVBaining']:
        a = generate_deni_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = False)
        a = generate_deni_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = True)

        a = generate_deni_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = False, resample_daf = 1)
        a = generate_deni_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = True, resample_daf = 1)
    for ascertainment_population in ['continent_europe', 'continent_sasia', 'continent_easia', 'subset_papuaExclFrancoisUVBaining']:
        a = generate_nean_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = False)
        a = generate_nean_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = True)

        a = generate_nean_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = False, resample_daf = 1)
        a = generate_nean_diagnostic(ascertainment_population = ascertainment_population, overwrite = False, africa_condition = True, resample_daf = 1)
