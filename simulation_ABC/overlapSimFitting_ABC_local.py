### ABC fitting of the introgression process. Runs locally but , and extended to include further calculations (haplotype lengths, population-wide overlaps)

import numpy as np
import subprocess
import tempfile
import sys
import os
import shutil
import copy
import time
import argparse

import msprime_functions_ABC

POPDATA_FOLDER = os.getcwd() + "/sample_lists/"
TMP_FOLDER = os.getcwd() + '/tmp_abc/'
OUT_FOLDER = os.getcwd() + '/out_abc'
MASK_FOLDER = os.getcwd() + '/Masks/'
BEDTOOLS_DIR = '' # Can e.g. point to specific version of bedtools or be left blank if bedtools is in a PATH directory


PIPE = ' | '
POP_DATA = {}
RUN_ID = 0


parser = argparse.ArgumentParser(description='Run an ABC sampling from priors for archaic introgression fitting on the Cambridge HPC cluster.')


parser.add_argument('--outfolder_base', metavar='-outfolder', dest='outfolder_base', type=str, nargs='?', default = OUT_FOLDER + '/IIM-T0=%d-T1=%d-a=%.3f-b=%.3f-c1=%.3f-c2=%.3f-N=%.3f-m1=%.6f-m2=%.6f/',
                    help='system location of the output folder root. The program will substitute in the relevant parameters and create this folder if required.')
parser.add_argument('--outfile_base', metavar='-outfile', dest='outfile_base', type=str, nargs='?', default = 'tmp%s_IIM-T0=%d-T1=%d-a=%.3f-b=%.3f-c1=%.3f-c2=%.3f-N=%.3f-m1=%.6f-m2=%.6f-ITx=%d-ITa=%d-ITb=%d-Ipx=%.6f-Ipa=%.6f-Ipb=%.6f-rep=%s',
                    help='system location of the output file root. The program will sub in the relevant parameters and write samples to this location.')

parser.add_argument('--T0_gens', '-T0', dest = 'T0_gens', type = int, action = 'store', nargs = '?', default = '4000',
                    help="time (generations, pastward) when the original population splits into two connected populations.")
parser.add_argument('--T1_gens', '-T1', dest = 'T1_gens', type = int, action = 'store', nargs = '?', default = '1000',
                    help="time (generations, pastward) when the migration between the two populations stops.")

parser.add_argument('--popsize_a', '-a', dest = 'popsize_a', type = float, action = 'store', nargs = '?', default = '10.0',
                    help="population size (in units relative to the size of the 'first' split population) of the original population.")
parser.add_argument('--popsize_b', '-b', dest = 'popsize_b', type = float, action = 'store', nargs = '?', default = '1.0',
                    help="population size (in units relative to the size of the 'first' split population) of the second split population.")
parser.add_argument('--popsize_c1', '-c1', dest = 'popsize_c1', type = float, action = 'store', nargs = '?', default = '2.0',
                    help="population size (in units relative to the size of the 'first' split population) of the first split population after migration stops.")
parser.add_argument('--popsize_c2', '-c2', dest = 'popsize_c2', type = float, action = 'store', nargs = '?', default = '2.0',
                    help="population size (in units relative to the size of the 'first' split population) of the second split population after migration stops.")
parser.add_argument('--popsize_N', '-N', dest = 'popsize_N', type = float, action = 'store', nargs = '?', default = '2000.0',
                    help="population size of the 'first' split population.")

parser.add_argument('--migration_1_AtoB', '-m1', dest = 'migration_1_AtoB', type = float, action = 'store', nargs = '?', default = '0.0001',
                    help="migration rate from the 'first' split population to the second, forward in time.")
parser.add_argument('--migration_2_BtoA', '-m2', dest = 'migration_2_BtoA', type = float, action = 'store', nargs = '?', default = '0.0001',
                    help="migration rate from the 'second' split population to the first, forward in time.")

parser.add_argument('--I_Tx_minmax', '-I_Tx_minmax', dest='I_Tx_minmax', type=int, action = 'store', nargs = '*', default = [4100,4100],
                    help="uniform prior over the initial introgression time (Ix) in generations pastward.")
parser.add_argument('--I_Ta_minmax', '-I_Ta_minmax', dest='I_Ta_minmax', type=int, action = 'store', nargs = '*', default = [1000,2000],
                    help="uniform prior over the time of the introgression into population A (Ia) in generations pastward.")
parser.add_argument('--I_Tb_minmax', '-I_Tb_minmax', dest='I_Tb_minmax', type=int, action = 'store', nargs = '*', default = [1000,2000],
                    help="uniform prior over the time of the introgression into population B (Ib) in generations pastward.")

parser.add_argument('--I_prop_x_minmax', '-I_prop_x_minmax', dest='I_prop_x_minmax', type=float, action = 'store', nargs = '*', default = [0.0,0.07],
                    help="by default, uniform prior over the initial introgression amount (Ix) in generations pastward.")
parser.add_argument('--I_prop_a_minmax', '-I_prop_a_minmax', dest='I_prop_a_minmax', type=float, action = 'store', nargs = '*', default = [0.0,0.07],
                    help="by default, uniform prior over the amount of the introgression into population A (Ia) in generations pastward.")
parser.add_argument('--I_prop_b_minmax', '-I_prop_b_minmax', dest='I_prop_b_minmax', type=float, action = 'store', nargs = '*', default = [0.0,0.04],
                    help="by default, uniform prior over the amount of the introgression into population B (Ib) in generations pastward.")

parser.add_argument('--I_prop_xab_minmax', '-I_prop_xab_minmax', dest='I_prop_xab_minmax', type=float, action = 'store', nargs = '*', default = [0.0,0.06],
                    help="usually ignored. By default, uniform prior over the total amount of the introgression experienced over events x, a and b; enabled using the samplemode switch.")

parser.add_argument('--samplemode', metavar='-smode', dest='samplemode', type=str, nargs='?', default = 'indep_intro_uniform',
                    help='mode of sampling. Often leave as default indep_intro_uniform which is just uniform sampling between min and max introgression pulse values. But options of total_intro which distributes a fixed total amount of introgression between the pulses, or indep_intro_expon_uniform which samples from an exponential.')
parser.add_argument('--nonuniform_intro', metavar='nonuni', dest = 'nonuniform_intro', type = str, nargs = '*', action = 'append',
                    help='when samplemode is indep_intro_expon_uniform (or, if implemented, other distributions) this can be used to specify which introgression events are not uniform.')

parser.add_argument('--iterations', metavar='-i', dest = 'iterations', type = int, action = 'store', nargs = '?', default = '1000',
                    help="how many iterations to run.")
parser.add_argument('--num_chromosomes', metavar='-chrom_n', dest = 'num_chromosomes', type = int, action = 'store', nargs = '?', default = '100',
                    help="number of chromosomes ie independent simulated chunks. Results (and speed) are a trade off between the number and size of chunks.")
parser.add_argument('--chrom_length', metavar='-chrom_L', dest = 'chrom_length', type = int, action = 'store', nargs = '?', default = '10000000',
                    help="length of chromosomes ie independent simulated chunks. Results (and speed) are a trade off between the number and size of chunks.")
parser.add_argument('--num_haps', metavar='-haps_n', dest = 'num_haps', type = int, action = 'store', nargs = '?', default = '50',
                    help="how many chromosome copies to simulate from each population. Results (and speed) are a trade off between the number and size of chunks.")
             
args = parser.parse_args()

#For testing

args.outfolder_base = os.getcwd() + '/out_abc/out_test/'
args.iterations = 10
args.num_chromosomes = 10
args.num_haps = 50
args.chrom_length = 1000000
args.I_Ta_minmax = [1700,1700]
args.I_Tb_minmax = [1700,1700]
args.samplemode = 'indep_intro_uniform'
args.I_prop_x_minmax = [0.0,0.03]
args.I_prop_a_minmax = [0.0,0.03]
args.I_prop_b_minmax = [0.0,0.03]


#args.samplemode = 'indep_intro_expon_uniform'
#args.nonuniform_intro = ['x', 'b']
#args.I_prop_x_minmax = [0.015,0.015]
#args.I_prop_a_minmax = [0.0,0.03]
#args.I_prop_b_minmax = [0.015,0.015]


def ABC_sample_IIM_multiAdmix(iterations = 100, outfolder_base = args.outfolder_base, outfile_base = 'tmp%s_IIM-T0=%d-T1=%d-a=%.3f-b=%.3f-c1=%.3f-c2=%.3f-N=%.3f-m1=%.6f-m2=%.6f-ITx=%d-ITa=%d-ITb=%d-Ipx=%.6f-Ipa=%.6f-Ipb=%.6f-rep=%s', T0_gens = 3502, T1_gens = 468, a = 7.37, b = 1.01, c1 = 47.29, c2 = 55.80, N = 2000.0, m1 = 0.00067, m2 = 0.00035, I_Tx_minmax = [3550,3550], I_Ta_minmax = [500,3000], I_Tb_minmax = [500,3000], I_tot_minmax = [0.0,0.06], I_prop_x_minmax = [0.0,0.06], I_prop_a_minmax = [0.0,0.06], I_prop_b_minmax = [0.0,0.06], I_T_resolution = None, I_prop_resolution = None, samplemode = 'total_intro', num_chromosomes = 10, chrom_length = 1e7, num_haps = 20):
    """
    Attempt an ABC fitting of the parameter values.
    I have a uniformly distributed amount of introgression (I_tot), which is then randomly split between the three events.
    """
    global RUN_ID
    for iteration in range(iterations):
        start_msprime = time.time()
        print("Iteration %d of %d at time %f" %(iteration, iterations, time.time()))
        I_Tx = (np.random.rand() * (I_Tx_minmax[1] - I_Tx_minmax[0])) + I_Tx_minmax[0]
        I_Ta = (np.random.rand() * (I_Ta_minmax[1] - I_Ta_minmax[0])) + I_Ta_minmax[0]
        I_Tb = (np.random.rand() * (I_Tb_minmax[1] - I_Tb_minmax[0])) + I_Tb_minmax[0]
        if samplemode == 'indep_intro_uniform':
            #This mode randomly and uniformly decides the introgression of each event independently
            I_prop_x = (np.random.rand() * (I_prop_x_minmax[1] - I_prop_x_minmax[0])) + I_prop_x_minmax[0]
            I_prop_a = (np.random.rand() * (I_prop_a_minmax[1] - I_prop_a_minmax[0])) + I_prop_a_minmax[0]
            I_prop_b = (np.random.rand() * (I_prop_b_minmax[1] - I_prop_b_minmax[0])) + I_prop_b_minmax[0]
        elif samplemode == 'total_intro':
            #This mode randomly and uniformly decides the total introgression.
            #It then randomly and fairly distributed this between the three introgression events.
            I_tot = (np.random.rand() * (I_tot_minmax[1] - I_tot_minmax[0])) + I_tot_minmax[0]
            I_tot_split = np.random.rand(3)
            I_tot_split = (I_tot_split / np.sum(I_tot_split)) * I_tot
            I_prop_x = I_tot_split[0]
            I_prop_a = I_tot_split[1]
            I_prop_b = I_tot_split[2]
        elif samplemode == 'indep_intro_expon_uniform':
            #This mode randomly decides the introgression of each event independently using a mix of exponential and uniform distributions.
            if ['x'] in args.nonuniform_intro:
                print("exponX")
                I_prop_x = np.random.exponential(I_prop_x_minmax[0])
            else:
                I_prop_x = (np.random.rand() * (I_prop_x_minmax[1] - I_prop_x_minmax[0])) + I_prop_x_minmax[0]
            if ['a'] in args.nonuniform_intro:
                print("exponA")
                I_prop_a = np.random.exponential(I_prop_a_minmax[0])
            else:
                I_prop_a = (np.random.rand() * (I_prop_a_minmax[1] - I_prop_a_minmax[0])) + I_prop_a_minmax[0]
            if ['b'] in args.nonuniform_intro:
                print("exponB")
                I_prop_b = np.random.exponential(I_prop_b_minmax[0])
            else:
                I_prop_b = (np.random.rand() * (I_prop_b_minmax[1] - I_prop_b_minmax[0])) + I_prop_b_minmax[0]
        #Discretise if that's the plan. Note that this may cause limited sampling at the boundaries if they are not chosen well.
        if type(I_T_resolution) != type(None):
            I_Tx = np.round(I_Tx / I_T_resolution) * I_T_resolution
            I_Ta = np.round(I_Ta / I_T_resolution) * I_T_resolution
            I_Tb = np.round(I_Tb / I_T_resolution) * I_T_resolution
        if type(I_prop_resolution) != type(None):
            I_prop_x = np.round(I_prop_x / I_prop_resolution) * I_prop_resolution
            I_prop_a = np.round(I_prop_a / I_prop_resolution) * I_prop_resolution
            I_prop_b = np.round(I_prop_b / I_prop_resolution) * I_prop_resolution

        try:
            outfolder = outfolder_base %(T0_gens, T1_gens, a, b, c1, c2, N, m1, m2)
        except TypeError:
            print("Unable to substitute run parameters into folder name, continuing.")
            outfolder = outfolder_base
        if os.path.exists(outfolder) == False:
            try:
                os.mkdir(outfolder)
            except FileExistsError:
                pass
        outfile = outfile_base %('%d', T0_gens, T1_gens, a, b, c1, c2, N, m1, m2, I_Tx, I_Ta, I_Tb, I_prop_x, I_prop_a, I_prop_b, "%s")
        rnd = np.random.randint(10000000000)
        while np.sum([os.path.exists(outfolder + outfile %(rnd, '0') + '.pop%s.blocks.source-%s.BED.gz' %(i, j)) for i, j in [['A', 'archA'], ['A', 'archB'], ['A', 'archX'], ['A','human'], ['B', 'archA'], ['B', 'archB'], ['B', 'archX'], ['B','human']]]) != 0:
            rnd = np.random.randint(10000000000)
        with open(outfolder + outfile %(rnd, '0') + '.popA.blocks.source-archA.BED.gz', 'wb') as f_tmp:
            pass
        RUN_ID = copy.copy(rnd)
        sim_outfile_name = outfolder + outfile %(rnd, '%d')
        msprime_functions_ABC.experiment_IIM_multipleintro(T0_gens = T0_gens,
                                     T1_gens = T1_gens,
                                     a = a,
                                     b = b,
                                     c1 = c1,
                                     c2 = c2,
                                     N = N,
                                     m1 = m1,
                                     m2 = m2,
                                     I_Tx = I_Tx,
                                     I_Ta = I_Ta,
                                     I_Tb = I_Tb,
                                     I_prop_x = I_prop_x,
                                     I_prop_a = I_prop_a,
                                     I_prop_b = I_prop_b,
                                     num_popA = num_haps,
                                     num_popB = num_haps,
                                     chromosomes = num_chromosomes,
                                     chrom_len = chrom_length,
                                     save_chroms = sim_outfile_name)

        inds = range(1, int(num_haps/2.0) + 1)
        for pop in ['popA', 'popB']:
            sim_infile_name = outfile %(rnd, '%s.' + pop)
            comb_outfile_name = outfile %(rnd, "COMBINED.%d.%d." + pop + '.blocks.source-%s.BED.gz')
            msprime_functions_ABC.combine_simulation_output(basefile_sims = sim_infile_name,
                                  basefolder_in = outfolder,
                                  basefolder_out = outfolder,
                                  reps_to_combine = range(num_chromosomes),
                                  ind_names = ['%s' %(i) for i in inds],
                                  vcf_reps_are_chroms = False,
                                  vcf_reps_are_chunks = False,
                                  vcf_reps_are_chunks_insets = None,
                                  rename_contig_and_filter_rep_vcfs = False,
                                  blocks = True,
                                  mismatch = False,
                                  mrca = False,
                                  sources = ['human', 'archA', 'archB', 'archX'])
            #Delete simulation files
            for rep in range(num_chromosomes):
                for file in [outfolder + sim_infile_name %(str(rep)) + '.blocks.source-%s.BED.gz' %(source) for source in ['human', 'archA', 'archB', 'archX']]:
                    os.remove(file)
            #Merge and delete excess files
            for ind in inds:
                for hap in range(1,3):
                    files_to_merge = [outfolder + comb_outfile_name %(ind, hap, source) for source in ['archA', 'archB', 'archX']]
                    msprime_functions_ABC.merge_beds(outfile = outfolder + comb_outfile_name[0:-3] %(ind, hap, ''.join(['archA', 'archB', 'archX'])), bed_list = files_to_merge, bedtools_dir = BEDTOOLS_DIR, tmp_folder = TMP_FOLDER)
                    for file in [outfolder + comb_outfile_name %(ind, hap, source) for source in ['human', 'archA', 'archB', 'archX']]:
                        os.remove(file)
            
        #Calculate the overlap
        overlap_file = outfile %(rnd, "COMBINED_pairwiseIndOver_simArchAll")
        introhapstats_file = outfile %(0, "COMBINED_pairwiseIndOver_simArchAll.introhapstats.csv")
        
        end_msprime = time.time()
        start_sumstat = time.time()

        #Work out which files to calculate pairwise overlap on
        
        POP_DATA['simpop_tmp'] = []
        POP_DATA['simpopA_tmp'] = []
        POP_DATA['simpopB_tmp'] = []
        popA_to_merge = []
        popB_to_merge = []
        tmp_file_list = []
        tmp_ind_list = []
        for pop in ['popA', 'popB']:
            for ind in inds:
                POP_DATA['simpop_tmp'].append(pop + '-%d' %(ind))
                for hap in [1,2]:
                    comb_file_name = outfile %(rnd, "COMBINED.%d.%d." + pop + '.blocks.source-archAarchBarchX.BED')
                    tmp_ind = pop + '-%d_%d' %(ind, hap)
                    tmp_file_list.append(outfolder + comb_file_name %(ind, hap))
                    tmp_ind_list.append(tmp_ind)
                    
                if pop == 'popA':
                    POP_DATA['simpopA_tmp'].append('popA-%d' %(ind))
                    popA_to_merge.append(outfolder + comb_file_name %(ind, 1))
                    popA_to_merge.append(outfolder + comb_file_name %(ind, 2))
                    
                elif pop == 'popB':
                    POP_DATA['simpopB_tmp'].append('popB-%d' %(ind))
                    popB_to_merge.append(outfolder + comb_file_name %(ind, 1))
                    popB_to_merge.append(outfolder + comb_file_name %(ind, 2))

        #Calculate pairwise overlap
        multipleAssess_pairwisePopOverlap_simpleABC(func_filelist = lambda : [np.array(tmp_file_list), np.array(tmp_ind_list)], outfile = outfolder + overlap_file, outfile_hapsumstats = outfolder + introhapstats_file, pops = ['simpop_tmp'], distance_approach = 'all', file_IDx_A = 0, genome_file_local = MASK_FOLDER + 'sim_%dchunksnum_%dMb.sizes' %(num_chromosomes, int(chrom_length/1000000)), genome_bed_local = MASK_FOLDER + 'sim_%dchunksnum_%dMb.sizes.bed' %(num_chromosomes, int(chrom_length/1000000)), genome_mask_local = MASK_FOLDER + 'sim_%dchunksnum_%dMb.mask' %(num_chromosomes, int(chrom_length/1000000)))
        
        #Merged population overlap
        msprime_functions_ABC.merge_beds_uncompressed(outfile = outfolder + outfile %(rnd, "COMBINED.merged.popA.blocks.source-archAarchBarchX.BED"), bed_list = popA_to_merge, bedtools_dir = BEDTOOLS_DIR, tmp_folder = TMP_FOLDER)
        msprime_functions_ABC.merge_beds_uncompressed(outfile = outfolder + outfile %(rnd, "COMBINED.merged.popB.blocks.source-archAarchBarchX.BED"), bed_list = popB_to_merge, bedtools_dir = BEDTOOLS_DIR, tmp_folder = TMP_FOLDER)
        popmerge_file_list = [outfolder + outfile %(rnd, "COMBINED.merged.popA.blocks.source-archAarchBarchX.BED"), outfolder + outfile %(rnd, "COMBINED.merged.popB.blocks.source-archAarchBarchX.BED")]
        popmerge_ind_list = ['simpopA_tmp--', 'simpopB_tmp--']
        POP_DATA['simpop_merged'] = ['simpopA_tmp', 'simpopB_tmp']

        multipleAssess_pairwisePopOverlap_simpleABC(func_filelist = lambda : [np.array(popmerge_file_list), np.array(popmerge_ind_list)], outfile = outfolder + outfile %(rnd, "COMBINED_pairwiseMergedPopOver_simArchAll"), outfile_hapsumstats = None, pops = ['simpop_merged'], distance_approach = 'all', file_IDx_A = 0, genome_file_local = MASK_FOLDER + 'sim_%dchunksnum_%dMb.sizes' %(num_chromosomes, int(chrom_length/1000000)), genome_bed_local = MASK_FOLDER + 'sim_%dchunksnum_%dMb.sizes.bed' %(num_chromosomes, int(chrom_length/1000000)), genome_mask_local = MASK_FOLDER + 'sim_%dchunksnum_%dMb.mask' %(num_chromosomes, int(chrom_length/1000000)))
        os.remove(outfolder + outfile %(rnd, "COMBINED.merged.popA.blocks.source-archAarchBarchX.BED"))
        os.remove(outfolder + outfile %(rnd, "COMBINED.merged.popB.blocks.source-archAarchBarchX.BED"))
        
        #Delete excess individual combined files
        for pop in ['popA', 'popB']:
            comb_file_name = outfile %(rnd, "COMBINED.%d.%d." + pop + '.blocks.source-archAarchBarchX.BED')
            for ind in inds:
                for hap in [1, 2]:
                    os.remove(outfolder + comb_file_name %(ind, hap))

        #Combine the individual overlap into populations and delete ind-vs-ind overlap files
        for approach in ['00','01','10','11']:
            pop_overlap_file = outfile %(0, "COMBINED_pairwiseIndOver_simArchAll") + '.%s.pops.%s.csv' %(approach, '%s')
            average_pairwise_ind_matrix_to_pop(in_matrix_ind = outfolder + overlap_file + '.%s.csv' %(approach),
                                   out_matrix_pop = outfolder + pop_overlap_file,
                                   pops = ['simpopA_tmp', 'simpopB_tmp'])
            os.remove(outfolder + overlap_file + '.%s.csv' %(approach))
        end_sumstat = time.time()
        
        #Combine the merged overlap into one file
        s11_s10_s01_s00 = []
        for approach in ['11', '10', '01', '00']:
            with open(outfolder + outfile %(rnd, "COMBINED_pairwiseMergedPopOver_simArchAll.%s.csv" %(approach)), 'rb') as f_in:
                headers = f_in.readline() #headers
                s11_s10_s01_s00.append(float(f_in.readline().decode().split(',')[2]))
            os.remove(outfolder + outfile %(rnd, "COMBINED_pairwiseMergedPopOver_simArchAll.%s.csv" %(approach)))
        popmerge_overlap_file = outfolder + outfile %(0, "COMBINED_pairwiseMergedPopOver_simArchAll") + '.s11s10s01s00.csv'
        with open(popmerge_overlap_file, 'wb') as f_out:
            f_out.write(','.join(['%.6f' %(i) for i in s11_s10_s01_s00]).encode())

        print("msprime step took %.2f seconds, sumstat step took %.1f seconds; total for one iteration %.1f seconds" %(end_msprime - start_msprime, end_sumstat - start_sumstat, end_sumstat - start_msprime))
        
    return None
    


def popAssess_pairwisePopOverlap_simpleFast(popbedtarg_a = [], popbedtarg_inds = [], popbedcomp_a = [], popbedcomp_inds = [], outfile_hapsumstats = None, pairwise_reciprocal = False, distance_approach = 'coverage', genome_file_local = None, genome_bed_local = None, genome_mask_local = None):
    """
    Assess the overlap between all individuals in one population and all individuals in a second population.
    If reciprocal is true then assess overlap both ways.
    
    This function is used to look at individual overlap; for example, if you want to calculate on all individuals you just pass all individuals in the popbedtarg_a and popbedcomp_a.

    This is a faster implementation to speed up for ABC.
    """
    
    #Create a local copy of the genome_files to avoid clashes when running the code on many cores
    genome_file_local = shutil.copy(src=genome_file_local, dst = TMP_FOLDER + '%d_' %(RUN_ID) + genome_file_local.split('/')[-1])
    genome_bed_local = shutil.copy(src=genome_bed_local, dst = TMP_FOLDER + '%d_' %(RUN_ID) + genome_bed_local.split('/')[-1])
    genome_mask_local = shutil.copy(src=genome_mask_local, dst = TMP_FOLDER + '%d_' %(RUN_ID) + genome_mask_local.split('/')[-1])

    popbedtarg_a = popbedtarg_a
    popbedcomp_a = popbedcomp_a
    
    ind_overlaps = {}
    pre_comp_time = time.time()
    
    #Calculate the sum of 1x for each individual. These can be used with 11 to calculate 01+10, which in turn lets us replace another pairwise comparison with a multi-comparison intersect, ie save a lot of time.
    overlap_command_1x = BEDTOOLS_DIR + 'bedtools intersect -a %s -b %s -g %s -wo' %(genome_bed_local, ' '.join(popbedtarg_a), genome_file_local) + PIPE + "awk '{sum1x[$4] += $8; count1x[$4] += 1}; $6 !~ /(0)/ && $7 !~ /(1000000)/ {sum1xNoEnd[$4] += $8; count1xNoEnd[$4] += 1} {sum1xNoEnd[$4] += 0; count1xNoEnd[$4] += 0}; $6 !~ /(0)/ && $7 !~ /(1000000)/ && $8 > 10000 {sum1xNoEndOver10k[$4] += $8; count1xNoEndOver10k[$4] += 1} {sum1xNoEndOver10k[$4] += 0; count1xNoEndOver10k[$4] += 0}; $6 !~ /(0)/ && $7 !~ /(1000000)/ && $8 > 20000 {sum1xNoEndOver20k[$4] += $8; count1xNoEndOver20k[$4] += 1} {sum1xNoEndOver20k[$4] += 0; count1xNoEndOver20k[$4] += 0}; END{ for (id in sum1x) { print id, sum1x[id], count1x[id], count1x[id] == 0 ? 0 : sum1x[id]/count1x[id], sum1xNoEnd[id], count1xNoEnd[id], count1xNoEnd[id] == 0 ? 0 : sum1xNoEnd[id]/count1xNoEnd[id], sum1xNoEndOver10k[id], count1xNoEndOver10k[id], count1xNoEndOver10k[id] == 0 ? 0 : sum1xNoEndOver10k[id]/count1xNoEndOver10k[id], sum1xNoEndOver20k[id], count1xNoEndOver20k[id], count1xNoEndOver20k[id] == 0 ? 0 : sum1xNoEndOver20k[id]/count1xNoEndOver20k[id]}}'" #"awk '{sum1x[$4] += $8}; END{ for (id in sum1x) { print id, sum1x[id]}}'"
    
    bedtools_process_1x = subprocess.Popen(args = overlap_command_1x, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    stdout, stderr = bedtools_process_1x.communicate(input = None)
    sum_1x = np.rot90(np.reshape((np.zeros(len(popbedtarg_a), dtype = int),np.arange(len(popbedtarg_a), dtype = int)), (2,len(popbedtarg_a))), k = 3)
    
    stats_1x = np.rot90(np.reshape((np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.zeros(len(popbedtarg_a), dtype = float),np.arange(len(popbedtarg_a), dtype = float)), (13,len(popbedtarg_a))), k = 3)
    for comp in np.array([line.split(b' ') for line in stdout.split(b'\n')][0:-1]):
        sum_1x[int(comp[0]) - 1][1] = int(comp[1])
        stats_1x[int(comp[0]) - 1][1:] = np.array(comp[1:], dtype = float)
    #Write out some haplotype length summary stats
    if type(outfile_hapsumstats) != type(None):
        #Write out some basic info on the chunk lengths
        with open(outfile_hapsumstats, 'wb') as f_sumintro_out:
            f_sumintro_out.write(b'hapID,sumIntro,countIntro,lengthIntro,sumIntroInternal,countIntroInternal,lengthIntroInternal,sumIntroInternalOver10k,countIntroInternalOver10k,lengthIntroInternalOver10k,sumIntroInternalOver20k,countIntroInternalOver20k,lengthIntroInternalOver20k\n')
            for line in stats_1x:
                f_sumintro_out.write(b'%d,%d,%d,%.1f,%d,%d,%.1f,%d,%d,%.1f,%d,%d,%.1f\n' %tuple(line))
        
    for ind_idx in range(len(popbedtarg_a)):
        ind = popbedtarg_inds[ind_idx]
        #print('Comparing %s with %s (excluding ideniticals)' %(ind, ' '.join(popbedcomp_inds)))

        #Calculating 11 overlaps
        num_comps = len(popbedcomp_a) - ind_idx - 1
        overlap_11 = np.rot90(np.reshape((np.zeros(num_comps, dtype = int),np.arange(num_comps, dtype = int)), (2,num_comps)), k = 3)
        if sum_1x[ind_idx][1] == 0:
            #No 11 overlaps possible as no 1x calls, leave everything as 0
            pass
        else:
            pre_comp_time = time.time()
            
            overlap_command_11 = BEDTOOLS_DIR + 'bedtools intersect -a %s -b %s -g %s -wao' %(popbedtarg_a[ind_idx], ' '.join(popbedcomp_a[i] for i in range(ind_idx + 1, len(popbedcomp_a))), genome_file_local) + PIPE + "awk '{sum11[$4] += $8}; END{ for (id in sum11) { print id, sum11[id]}}'"
            if len([popbedcomp_a[i] for i in range(ind_idx + 1, len(popbedcomp_a))]) == 1:
                overlap_command_11 = BEDTOOLS_DIR + 'bedtools intersect -a %s -b %s -g %s -wao' %(popbedtarg_a[ind_idx], ' '.join(popbedcomp_a[i] for i in range(ind_idx + 1, len(popbedcomp_a))), genome_file_local) + PIPE + "awk '{sum11 += $7}; END{ print 1, sum11}'"
            bedtools_process_11 = subprocess.Popen(args = overlap_command_11, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
            stdout, stderr = bedtools_process_11.communicate(input = None)
            for comp in np.array([line.split(b' ') for line in stdout.split(b'\n')][0:-1]):
                if comp[0] != b'.':
                    overlap_11[int(comp[0]) - 1][1] = int(comp[1])
        
        for comp_idx in range(ind_idx + 1, len(popbedcomp_a)):
            comp_ind = popbedcomp_inds[comp_idx]
            #print(ind, comp_ind)

            # Compare the two individuals based on coverage or jaccard distance.
            if distance_approach == 'coverage':
                #Coverage is 11 / [11 + 10] and is asymmetrical
                overlap_command = BEDTOOLS_DIR + 'bedtools coverage -a %s -b %s -g %s' %(popbedtarg_a[ind_idx], popbedcomp_a[comp_idx], genome_file_local) + PIPE + "awk '{SUMcover+=$5;SUMtot+=$6}END{print SUMcover/SUMtot}'"
                output_transform = lambda a : np.array([float(a)])
            elif distance_approach == 'jaccard':
                #Jaccard is 11 / [11 + 10 + 01] and so is symmetrical
                overlap_command = BEDTOOLS_DIR + 'bedtools jaccard -a %s -b %s -g %s' %(popbedtarg_a[ind_idx], popbedcomp_a[comp_idx], genome_file_local)
                output_transform = lambda a : np.array([float(a.split('\n')[1].split('\t')[2])])
            elif distance_approach == 'all':
                
                result_11_10_01_00 = []
                result_11_10_01_00.append(overlap_11[comp_idx - ind_idx - 1][1])
                result_11_10_01_00.append(sum_1x[ind_idx][1] - overlap_11[comp_idx - ind_idx - 1][1]) #1x - 11 ie calls in ind_idx that are not also in comp_idx
                result_11_10_01_00.append(sum_1x[comp_idx][1] - overlap_11[comp_idx - ind_idx - 1][1]) #x1 - 11 ie calls in comp_idx that are not also in ind_idx
                result_11_10_01_00.append((args.num_chromosomes * args.chrom_length) - np.sum(result_11_10_01_00)) #Ie everything that isn't 11, 10 or 01
                result_11_10_01_00 = np.array(result_11_10_01_00)
                
                try:
                    ind_overlaps[ind][comp_ind] = result_11_10_01_00
                except KeyError:
                    ind_overlaps[ind] = {}
                    ind_overlaps[ind][comp_ind] = result_11_10_01_00
            
            if distance_approach != 'all':    
                bedtools_process_coverage = subprocess.Popen(args = overlap_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                stdout, stderr = bedtools_process_coverage.communicate(input = None)
                if stderr is not '':
                    print(overlap_command)
                    print(stdout)
                    raise RuntimeError(stderr)
                try:
                    ind_overlaps[ind][comp_ind] = output_transform(stdout)
                except KeyError:
                    ind_overlaps[ind] = {}
                    ind_overlaps[ind][comp_ind] = output_transform(stdout)
            if pairwise_reciprocal == True:
                if distance_approach == 'coverage':
                    overlap_command = BEDTOOLS_DIR + 'bedtools coverage -a %s -b %s -g %s' %(popbedcomp_a[comp_idx], popbedtarg_a[ind_idx], genome_file_local) + PIPE + "awk '{SUMcover+=$5;SUMtot+=$6}END{print SUMcover/SUMtot}'"
                elif distance_approach == 'jaccard':
                    # The actual values of the reciprocal should be idenitical. They are very close, but not quite identical.
                    # I therefore include both, despite the speed cost.
                    overlap_command = BEDTOOLS_DIR + 'bedtools jaccard -a %s -b %s -g %s' %(popbedcomp_a[comp_idx], popbedtarg_a[ind_idx], genome_file_local)
                elif distance_approach == 'all':
                    # The overlaps are symmetrical for the reciprocal 'all' case! So no need to calculate again.
                    result_11_10_01_00_reciprocal = np.array([ind_overlaps[ind][comp_ind][0], ind_overlaps[ind][comp_ind][2], ind_overlaps[ind][comp_ind][1], ind_overlaps[ind][comp_ind][3]])
                    try:
                        ind_overlaps[comp_ind][ind] = result_11_10_01_00_reciprocal
                    except KeyError:
                        ind_overlaps[comp_ind] = {}
                        ind_overlaps[comp_ind][ind] = result_11_10_01_00_reciprocal

                if distance_approach != 'all':
                    bedtools_process_coverage = subprocess.Popen(args = overlap_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                    stdout, stderr = bedtools_process_coverage.communicate(input = None)
                    if stderr is not '':
                        print(overlap_command)
                        print(stdout)
                        raise RuntimeError(stderr)
                    try:
                        ind_overlaps[comp_ind][ind] = output_transform(stdout)
                    except KeyError:
                        ind_overlaps[comp_ind] = {}
                        ind_overlaps[comp_ind][ind] = output_transform(stdout)
            # print ind, comp_ind, ind_overlaps[comp_ind][ind]
            # print comp_ind, ind, ind_overlaps[ind][comp_ind]
        #print ind_overlaps[ind]

    os.remove(genome_file_local)
    os.remove(genome_bed_local)
    os.remove(genome_mask_local)
    
    return ind_overlaps


def multipleAssess_pairwisePopOverlap_simpleABC(func_filelist, outfile, outfile_hapsumstats = None, pops = ['sims_all75'], distance_approach = 'coverage', file_IDx_A = 0, genome_file_local = None, genome_bed_local = None, genome_mask_local = None):
    """
    ADAPTED FROM bedtools_archaic_v1
    
    Ask for overlap between sets of individuals using popAssess_pairwisePopOverlap_simpleFast.
    
    In this case, I also want to average over _1 and _2 of individuals.
    """
    pairwise_starttime = time.time()
    assess_files, ind_names = func_filelist()

    #if average_approach is 'population':
    #    print(RuntimeWarning("NOTE that population mode involves comparing all individuals in each included population with each included population, and then averaging, not including self-comparisons. In this case a list of population is appropriate."))
    #    if len(pops) <= 1:
    #        print("More than one pop needed; try pops = ['pop_Sulaw_K', 'pop_Sulaw_M', 'subset_lembAll', 'pop_Fl_Ram', 'pop_Fl_Ben', 'pop_Fl_Cib', 'pop_Tanim', 'pop_Alor', 'pop_Kei'] + ['pop_Pap_MPI', 'pop_Pap_Bai', 'selection_Papua_Highland_1', 'selection_NG2', 'selection_Korowai_and_East_Sepik', 'subset_papuaMalaspinas_pop1tmp', 'subset_papuaMalaspinas_pop2tmp', 'subset_papuaMalaspinas_pop3tmp', 'subset_papuaMalaspinas_pop4tmp', 'subset_papuaMalaspinas_pop5tmp'] ?")
    #        raise RuntimeError()
    #elif average_approach is 'individual':
    #    print(RuntimeWarning("NOTE that individual mode involves comparing all individuals in the population with all other individuals. A complete individual list is used. Individual populations are _1 vs _2 and _2 vs _1 (not self-comparison)."))
    order = np.sort(np.array(pops))
    all_results = {}
    av_results = {}
    for pop_targ_idx in np.arange(len(order)):
        #go through each
        pop_targ = order[pop_targ_idx]
        #print("Calculating overlaps for target population %s" %(pop_targ))
        pop_results = {}
        poptarget_assess = []
        poptarget_inds = []
        for ind in range(len(ind_names)):
            if ind_names[ind][0:-2] in POP_DATA[pop_targ]:
                poptarget_assess.append(assess_files[ind])
                poptarget_inds.append(ind_names[ind])
        poptarget_assess = np.array(poptarget_assess)
        poptarget_inds = np.array(poptarget_inds)
        #Ultimately, I want:
        #IndA sharingPopA sharingPopB sharingPopC ...
        for pop_comp_idx in np.arange(len(order) - pop_targ_idx) + pop_targ_idx:
            pop_comp = order[pop_comp_idx]
            print("Comparing with comparison population %s" %(pop_comp))
            popcomp_assess = []
            popcomp_inds = []
            for ind in range(len(ind_names)):
                if ind_names[ind][0:-2] in POP_DATA[pop_comp]:
                    popcomp_assess.append(assess_files[ind])
                    popcomp_inds.append(ind_names[ind])
            popcomp_assess = np.array(popcomp_assess)
            popcomp_inds = np.array(popcomp_inds)
            #popcomp_assess = popcomp_assess[0:6]
            #popcomp_inds = popcomp_inds[0:6]
            #if average_approach == 'population':
            #    pairwise_reciprocal = True if pop_comp != pop_targ else False
            #elif average_approach == 'individual':
            #    pairwise_reciprocal = True
            #else:
            #    raise RuntimeError()
            #print(average_approach)
            pop_results_targ = popAssess_pairwisePopOverlap_simpleFast(popbedtarg_a = poptarget_assess,#[::,0],
                                                    #popbedtarg_b = poptarget_assess[::,1],
                                                    popbedtarg_inds = poptarget_inds,
                                                    popbedcomp_a = popcomp_assess,#[::,0],
                                                    #popbedcomp_b = popcomp_assess[::,1],
                                                    popbedcomp_inds = popcomp_inds,
                                                    outfile_hapsumstats = outfile_hapsumstats,
                                                    pairwise_reciprocal = True,
                                                    distance_approach = distance_approach,
                                                    #func_conversion = func_conversion,
                                                    #resamples = 100,
                                                    #genome_mask = None,
                                                    genome_file_local = genome_file_local,
                                                    genome_bed_local = genome_bed_local,
                                                    genome_mask_local = genome_mask_local)
            for ind in pop_results_targ.keys():
                if ind[0:-2] not in all_results.keys():
                    all_results[ind[0:-2]] = {}
                for comp in pop_results_targ[ind].keys():
                    if comp[0:-2] not in all_results[ind[0:-2]].keys():
                        all_results[ind[0:-2]][comp[0:-2]] = [pop_results_targ[ind][comp]]
                    else:
                        all_results[ind[0:-2]][comp[0:-2]].append(pop_results_targ[ind][comp])
            #print all_results
    print("COMPLETED! Writing.")
    f_list = [outfile] if distance_approach != 'all' else [outfile + '.11.csv', outfile + '.10.csv', outfile + '.01.csv', outfile + '.00.csv']
    for f_out_idx in range(len(f_list)):
        f_out = open(f_list[f_out_idx], 'wb')
        av_results = {}
        for ind in all_results.keys():
            av_results[ind] = {}
            for comp in all_results[ind].keys():
                av_results[ind][comp] = np.average(np.array(all_results[ind][comp])[::,f_out_idx])
        ind_order = np.sort(list(all_results.keys()))
        #print(ind_order)
        f_out.write((','.join(['IND'] + list(ind_order)) + '\n').encode())
        for ind in ind_order:
            line_to_write = [ind]
            for comp_ind in ind_order:
                try:
                    line_to_write.append('%.3f' %(av_results[ind][comp_ind]))
                except:
                    line_to_write.append('nan')
            f_out.write((','.join(line_to_write) + '\n').encode())
        f_out.close()
    print("Total pairwise comparison time", time.time() - pairwise_starttime)
    return all_results


def average_pairwise_ind_matrix_to_pop(in_matrix_ind = '/home/guy/Dropbox/Papers/worldwideDeniNeanSharing/IIM_overlap/data/pairwiseIndOver_coverage_allNonSSAf_deniHCSS35unique.11.csv', out_matrix_pop = '/home/guy/Dropbox/Papers/worldwideDeniNeanSharing/IIM_overlap/data/pairwiseIndOver_coverage_allNonSSAf_deniHCSS35unique.11.pops.%s.csv', pops = ['continent_america', 'continent_e_isea', 'continent_easia', 'continent_europe', 'continent_sasia', 'continent_seasia', 'continent_siberia', 'continent_w_isea', 'subset_papuaExclFrancois']):
    """
    From overlapTheoryFitting.py
    
    Read in a matix of individual pairwise comparisons and write out the mean and the standard deviation of
    pairwise comparisons between a given set of populations.
    """
    matrix_ind = {}
    with open(in_matrix_ind, 'rb') as f_in:
        headers = f_in.readline()[0:-1].split(b',')
        for line in f_in:
            split_line = line[0:-1].split(b',')
            ind = split_line[0]
            for compind_idx in range(1, len(split_line)):
                try:
                    matrix_ind[ind][headers[compind_idx]] = float(split_line[compind_idx])
                except KeyError:
                    matrix_ind[ind] = {}
                    matrix_ind[ind][headers[compind_idx]] = float(split_line[compind_idx])
        matrix_pop_av = {}
        matrix_pop_std = {}
        for pop_targ in pops:
            for pop_comp in pops:
                #print(pop_targ, pop_comp)
                inds_targ = POP_DATA[pop_targ]
                inds_comp = POP_DATA[pop_comp]
                to_compare = []
                for ind_targ in inds_targ:
                    for ind_comp in inds_comp:
                        to_compare.append(matrix_ind[ind_targ.encode()][ind_comp.encode()])
                        try:
                            matrix_pop_av[pop_targ][pop_comp] = np.average(to_compare)
                        except KeyError:
                            matrix_pop_av[pop_targ] = {}
                            matrix_pop_av[pop_targ][pop_comp] = np.average(to_compare)
                        try:
                            matrix_pop_std[pop_targ][pop_comp] = np.std(to_compare)
                        except KeyError:
                            matrix_pop_std[pop_targ] = {}
                            matrix_pop_std[pop_targ][pop_comp] = np.std(to_compare)
        with open(out_matrix_pop %('av'), 'wb') as f_out:
            headers = ["POP"] + pops
            f_out.write((','.join(headers) + '\n').encode())
            for pop_targ in pops:
                line_to_write = [pop_targ]
                for pop_comp in pops:
                    line_to_write.append('%.3f' %(matrix_pop_av[pop_targ][pop_comp]))
                f_out.write((','.join(line_to_write) + '\n').encode())
        with open(out_matrix_pop %('std'), 'wb') as f_out:
            headers = ["POP"] + pops
            f_out.write((','.join(headers) + '\n').encode())
            for pop_targ in pops:
                line_to_write = [pop_targ]
                for pop_comp in pops:
                    line_to_write.append('%.3f' %(matrix_pop_std[pop_targ][pop_comp]))
                f_out.write((','.join(line_to_write) + '\n').encode())
    return None

if True:
    #Create the out and tmp directories if they don't exist
    if os.path.exists(TMP_FOLDER) == False:
        try:
            os.mkdir(TMP_FOLDER)
        except FileExistsError:
            pass
    if os.path.exists(OUT_FOLDER) == False:
        try:
            os.mkdir(OUT_FOLDER)
        except FileExistsError:
            pass
    a = ABC_sample_IIM_multiAdmix(iterations = args.iterations, outfolder_base = args.outfolder_base, outfile_base = args.outfile_base, T0_gens = args.T0_gens, T1_gens = args.T1_gens, a = args.popsize_a, b = args.popsize_b, c1 = args.popsize_c1, c2 = args.popsize_c2, N = args.popsize_N, m1 = args.migration_1_AtoB, m2 = args.migration_2_BtoA, I_Tx_minmax = args.I_Tx_minmax, I_Ta_minmax = args.I_Ta_minmax, I_Tb_minmax = args.I_Tb_minmax, I_prop_x_minmax = args.I_prop_x_minmax, I_prop_a_minmax = args.I_prop_a_minmax, I_prop_b_minmax = args.I_prop_b_minmax, I_tot_minmax = args.I_prop_xab_minmax, samplemode = args.samplemode, num_chromosomes = args.num_chromosomes, chrom_length = args.chrom_length, num_haps = args.num_haps)
