## This is code for processing the simulation output.
## It is run using a directory of simulations as the input. It creates a CSV
## that is easily used for later ABC fitting.

## GSJ 14-03-23, modified from original code to accept command line inputs and
## to work with the StdDev data too.
## GSJ 13-04-23, modified to work with population-merged overlap files.

import numpy as np
import os
import argparse


parser = argparse.ArgumentParser(description='Process the ABC simulation outputs into an easy-to-work with CSV file for fitting.')

parser.add_argument('--simfolder', metavar='-s_in', dest='simfolder', type=str, nargs='?', default = os.getcwd() + '/out_abc/out_test/',
                    help='system location of the folder containing the ABC simulation outputs. This will contain a lot of files ending in .pops.av.csv and .pops.std.csv ')
parser.add_argument('--obsfolder', metavar='-o_in', dest='obsfolder', type=str, nargs='?', default = '/rds/user/gsj22/rds-mobile-Gwxt74iudAM/Archaic_Sharing/ABC/Obs_rawfiles/',
                    help='system location of the folder containing the observed data to fit against. This is used to calculate the distance metric.')
parser.add_argument('--outfile_base', metavar='-out', dest='outfile_base', type=str, nargs='?', default = os.getcwd() + '/grid_abc/out_test',
                    help='system location of the output files. Note that intermediate files are stored, leave off the .csv .')

parser.add_argument('--calc_dist', metavar='-dist', dest='calc_dist', type=bool, nargs='?', default = False,
                    help='should the program calculate a distance metric? If not, I obsfolder and populations are not used.')
parser.add_argument('--grid_mode', metavar='-mode', dest='grid_mode', type=str, action = 'store', nargs = '?', default = 'av_only',
                    help="calculate the grid based on average haplotype sharing only ('avonly') or based on all data recorded from the sims ('all'). Note that power corrections are not applied to all outputs, and some outputs may be difficult to compare between simulated and observed data.")


parser.add_argument('--archmode', '-arch', dest = 'archmode', type = str, action = 'store', nargs = '?', default = 'deni',
                    help="run the program in 'deni' (Denisovan only) or 'nean' (Neanderthal only) mode.")
parser.add_argument('--arch_power', '-pow', dest='arch_power', type=float, action = 'store', nargs = '?', default = '0.29',
                    help="assumed power to detect haplotypes using the original methods applied. Essentially, [observed total introgression based on haplotype method/ assumed real total introgression]. Leave as default for deni analysis based on GSJ haplotypes, set to 0.39 for Nean analysis based on GSJ haplotypes.")
parser.add_argument('--arch_corr', '-cor', dest='arch_corr', type=float, action = 'store', nargs = '*', default = [0.0,1.0,3],
                    help="assumed correlation in power to identify archaic haplotypes between samples. 0 means random, 1 means perfect correlation. Leave as default to explore this assumption.")

parser.add_argument('--populations', '-pops', dest='populations', type=str, action = 'store', nargs = '*', default = ['continent_europe','continent_easia'],
                    help="populations that we are fitting KL divergence of the statistics to. Note this shouldn't be used if we are just summarising the stats for later processing, eg fitting on derived stats.")

parser.add_argument('--derived_parameter_targs', '--d_targ', dest='derived_parameter_targs', type = str, action = 'append', nargs = '*', default = [],
                    help="append a list of targets to our derived parameter calculation operation. The default tells the program to apply the first derived_parameter_action operation on the three I_prop values.")
parser.add_argument('--derived_parameter_operations', '--d_op', dest='derived_parameter_operations', type = str, action = 'append', nargs = '*', default = [],
                    help="append a list of operations to apply to the derived parameter targets. Default operations are sum (on all inputs), divide (first input by second input), multiply (all inputs) and subtract (first input minus second input). The function also accepts lambda functions, which operate on an array eg lambda a : np.product(a) for multiply. Make sure the expected number of input parameters are given in each case.")
parser.add_argument('--derived_parameter_names', '--d_names', dest='derived_parameter_names', type = str, action = 'append', nargs = '*', default = [],
                    help="append a derived parameter name to our list of derived parameters. Each new parameter requires a name.")


args = parser.parse_args()

##Testing
"""
args.derived_parameter_targs = [['I_prop_x','I_prop_a','I_prop_b'],
                                ['s11_AA', 's11_AA', 's10_AA'],
                                ['s11_BB', 's11_BB', 's10_BB'],
                                ['s11_AB', 's11_AB', 's10_AB'],
                                ['s11_AB', 's11_AB', 's01_AB'], #Note using AB as BA not saved and symmetrical
                                ['s11_AA', 's11_AA', 's10_AA', 's01_AA'],
                                ['s11_BB', 's11_BB', 's10_BB', 's01_BB'],
                                ['s11_AB', 's11_AB', 's10_AB', 's01_AB'], #Symmetrical to BA and BA not saved so use AB.
                                ['coverage_AB', 'coverage_BA']
                                ]
args.derived_parameter_operations = [['sum'],
                                     ['lambda x : x[0] / float(x[1] + x[2])'],
                                     ['lambda x : x[0] / float(x[1] + x[2])'],
                                     ['lambda x : x[0] / float(x[1] + x[2])'],
                                     ['lambda x : x[0] / float(x[1] + x[2])'],
                                     ['lambda x : x[0] / float(x[1] + x[2] + x[3])'],
                                     ['lambda x : x[0] / float(x[1] + x[2] + x[3])'],
                                     ['lambda x : x[0] / float(x[1] + x[2] + x[3])'],
                                     ['lambda x : x[0] - x[1]'],
                                     ]
args.derived_parameter_names = [['I_prop_total'],
                                ['coverage_AA'],
                                ['coverage_BB'],
                                ['coverage_AB'],
                                ['coverage_BA'],
                                ['jaccard_AA'],
                                ['jaccard_BB'],
                                ['jaccard_AB_BA'],
                                ['imbalance_metric']
                                ]
"""

def create_lambda_with_globals(s):
    return eval(s, globals())

##Test generally. Decide on ratio stats etc. Try running...

def construct_grid_from_ABCdir_avOnly(ABCdir = os.getcwd() + '/out_abc/out_test/', gridfile_out = '', power_min_max_res = [1.0,1.0,1], powerCorr_min_max_res = [0.0,0.0,1]):
    power_grid = np.linspace(power_min_max_res[0],power_min_max_res[1],int(power_min_max_res[2]))
    powerCorr_grid = np.linspace(powerCorr_min_max_res[0],powerCorr_min_max_res[1],int(powerCorr_min_max_res[2]))
    
    sim_filebase_set = set()
    sim_files = os.listdir(ABCdir)
    sim_files_set = set(sim_files)
    for file in sim_files:
        skip = False
        if file[0:5] != 'tmp0_':
            pass
        elif 'av.csv' in file:
            filetest_start = file.split('COMBINED_pairwiseIndOver_simArchAll.')[0] + 'COMBINED_pairwiseIndOver_simArchAll.'
            filetest_end = '.'.join(file.split('COMBINED_pairwiseIndOver_simArchAll.')[1].split('.')[1:])
            filetest_name = filetest_start + '%s.' + filetest_end
            #print(filetest_name)
            if filetest_name not in sim_filebase_set:
                #print([ABCdir + filetest_name %(i) for i in ['11', '00', '10', '01']])
                if np.sum([filetest_name %(i) in sim_files_set for i in ['11', '00', '10', '01']]) == 4:
                    sim_filebase_set.add(ABCdir + filetest_name)
        else:
            pass
    print("Read in %d states" %(len(sim_filebase_set)))
    sim_results_dict = {}
    for filetest_name in list(sim_filebase_set):
        #print(filetest_name)
        #Work out what the parameters are and record the output
        ITx = filetest_name.split('ITx=')[1].split('-')[0]
        ITa = filetest_name.split('ITa=')[1].split('-')[0]
        ITb = filetest_name.split('ITb=')[1].split('-')[0]
        Ipx = filetest_name.split('Ipx=')[1].split('-')[0]
        Ipa = filetest_name.split('Ipa=')[1].split('-')[0]
        Ipb = filetest_name.split('Ipb=')[1].split('-')[0]
        s11_s10_s01_s00_AABBAB = []
        for approach in ['11','10','01','00']:
            over_approach = []
            with open(filetest_name %(approach)) as f_popover:
                for line in f_popover:
                    split_line = line[0:-1].split(',')
                    if split_line[0] == 'simpopA_tmp':
                        over_approach.append(float(split_line[1])) #AA
                        over_approach.append(float(split_line[2])) #AB
                    elif split_line[0] == 'simpopB_tmp':
                        over_approach.append(float(split_line[2])) #BB
            s11_s10_s01_s00_AABBAB.append(over_approach)
        s11_s10_s01_s00 = np.array([[s11_s10_s01_s00_AABBAB[i][j] for i in [0,1,2,3]] for j in [0,2,1]])
        s11_s10_s01_s00 = s11_s10_s01_s00 / np.vstack(np.sum(s11_s10_s01_s00, 1)) #Proportions?
        #print(s11_s10_s01_s00)
        for power in power_grid:
            for powerCorr in powerCorr_grid:
                sim_results_dict[','.join([ITx, Ipx, ITa, Ipa, ITb, Ipb, '%.3f' %(power), '%.3f' %(powerCorr)])] = [modify_11_10_01_00_on_binom(popSamplePair, p = 1.0 - power, q = 1.0 - power, rho = powerCorr) for popSamplePair in s11_s10_s01_s00]
    print("Total %d states, writing" %(len(sim_results_dict)))
    with open(gridfile_out, 'wb') as f_out:
        headers = ['I_Tx', 'I_prop_x', 'I_Ta', 'I_prop_a', 'I_Tb', 'I_prop_b', 'power', 'powerCorrRho', 's11_AA', 's10_AA', 's01_AA', 's00_AA', 's11_BB', 's10_BB', 's01_BB', 's00_BB', 's11_AB', 's10_AB', 's01_AB', 's00_AB']
        f_out.write((','.join(headers) + '\n').encode())
        for key in np.sort(list(sim_results_dict.keys())):
            line_to_write = ','.join([key, ','.join([','.join(['%.6f' %(j) for j in i ]) for i in sim_results_dict[key]])])
            f_out.write((line_to_write + '\n').encode())
    return len(sim_filebase_set)

def construct_grid_from_ABCdir_allData(ABCdir = os.getcwd() + '/out_abc/out_test/', gridfile_out = '', power_min_max_res = [1.0,1.0,1], powerCorr_min_max_res = [0.0,0.0,1]):
    ## This is 'in development' - I'm not sure what the impact of power would be on these results, so I'd be very wary of using
    ## them in the fittings.
    power_grid = np.linspace(power_min_max_res[0],power_min_max_res[1],int(power_min_max_res[2]))
    powerCorr_grid = np.linspace(powerCorr_min_max_res[0],powerCorr_min_max_res[1],int(powerCorr_min_max_res[2]))
    
    sim_filebase_set_av = set()
    sim_filebase_set_std = set()
    sim_filebase_set_haps = set()
    sim_filebase_set_merged = set()
    
    sim_files = os.listdir(ABCdir)
    sim_files_set = set(sim_files)
    for file in sim_files:
        skip = False
        if file[0:5] != 'tmp0_':
            pass
        elif 'std.csv' in file:
            filetest_start = file.split('COMBINED_pairwiseIndOver_simArchAll.')[0] + 'COMBINED_pairwiseIndOver_simArchAll.'
            filetest_end = '.'.join(file.split('COMBINED_pairwiseIndOver_simArchAll.')[1].split('.')[1:])
            filetest_name = filetest_start + '%s.' + filetest_end
            if filetest_name not in sim_filebase_set_std:
                if np.sum([filetest_name %(i) in sim_files_set for i in ['11', '00', '10', '01']]) == 4:
                    sim_filebase_set_std.add(ABCdir + filetest_name)
        elif 'av.csv' in file:
            filetest_start = file.split('COMBINED_pairwiseIndOver_simArchAll.')[0] + 'COMBINED_pairwiseIndOver_simArchAll.'
            filetest_end = '.'.join(file.split('COMBINED_pairwiseIndOver_simArchAll.')[1].split('.')[1:])
            filetest_name = filetest_start + '%s.' + filetest_end
            #print(filetest_name)
            if filetest_name not in sim_filebase_set_av:
                if np.sum([filetest_name %(i) in sim_files for i in ['11', '00', '10', '01']]) == 4:
                    sim_filebase_set_av.add(ABCdir + filetest_name)
        elif 'introhapstats.csv' in file:
            if file not in sim_filebase_set_haps:
                sim_filebase_set_haps.add(ABCdir + file)
        elif 'pairwiseMergedPopOver_simArchAll.s11s10s01s00.csv' in file:
            #Population-merged overlap files
            if file not in sim_filebase_set_merged:
                sim_filebase_set_merged.add(ABCdir + file)
        else:
            pass
    ###There is going to be a problem here if I am missing some files for some sims.
    print("Read in %d states" %(len(sim_filebase_set_av)))
    sim_results_dict_av = {}
    sim_results_dict_std = {}
    sim_results_dict_haps = {}
    sim_results_dict_merged = {}
    
    total_bp = None
    for filetest_name in list(sim_filebase_set_av):
        #print(filetest_name)
        #Work out what the parameters are and record the output
        ITx = filetest_name.split('ITx=')[1].split('-')[0]
        ITa = filetest_name.split('ITa=')[1].split('-')[0]
        ITb = filetest_name.split('ITb=')[1].split('-')[0]
        Ipx = filetest_name.split('Ipx=')[1].split('-')[0]
        Ipa = filetest_name.split('Ipa=')[1].split('-')[0]
        Ipb = filetest_name.split('Ipb=')[1].split('-')[0]
        s11_s10_s01_s00_AABBAB = []
        for approach in ['11','10','01','00']:
            over_approach = []
            with open(filetest_name %(approach)) as f_popover:
                for line in f_popover:
                    split_line = line[0:-1].split(',')
                    if split_line[0] == 'simpopA_tmp':
                        over_approach.append(float(split_line[1])) #AA
                        over_approach.append(float(split_line[2])) #AB
                    elif split_line[0] == 'simpopB_tmp':
                        over_approach.append(float(split_line[2])) #BB
            s11_s10_s01_s00_AABBAB.append(over_approach)
        s11_s10_s01_s00 = np.array([[s11_s10_s01_s00_AABBAB[i][j] for i in [0,1,2,3]] for j in [0,2,1]])

        if type(total_bp) == type(None):
            total_bp = np.sum(s11_s10_s01_s00, 1)[0]
        else:
            assert total_bp == np.sum(s11_s10_s01_s00, 1)[0]

        assert np.sum(np.vstack(np.sum(s11_s10_s01_s00, 1))) == (3 * total_bp)
        
        s11_s10_s01_s00 = s11_s10_s01_s00 / np.vstack(np.sum(s11_s10_s01_s00, 1))
        #print(s11_s10_s01_s00)
        for power in power_grid:
            for powerCorr in powerCorr_grid:
                sim_results_dict_av[','.join([ITx, Ipx, ITa, Ipa, ITb, Ipb, '%.3f' %(power), '%.3f' %(powerCorr)])] = [modify_11_10_01_00_on_binom(popSamplePair, p = 1.0 - power, q = 1.0 - power, rho = powerCorr) for popSamplePair in s11_s10_s01_s00]

    for filetest_name in list(sim_filebase_set_std):
        #print(filetest_name)
        #Work out what the parameters are and record the output
        ITx = filetest_name.split('ITx=')[1].split('-')[0]
        ITa = filetest_name.split('ITa=')[1].split('-')[0]
        ITb = filetest_name.split('ITb=')[1].split('-')[0]
        Ipx = filetest_name.split('Ipx=')[1].split('-')[0]
        Ipa = filetest_name.split('Ipa=')[1].split('-')[0]
        Ipb = filetest_name.split('Ipb=')[1].split('-')[0]
        s11_s10_s01_s00_AABBAB = []
        for approach in ['11','10','01','00']:
            over_approach = []
            with open(filetest_name %(approach)) as f_popover:
                for line in f_popover:
                    split_line = line[0:-1].split(',')
                    if split_line[0] == 'simpopA_tmp':
                        over_approach.append(float(split_line[1])) #AA
                        over_approach.append(float(split_line[2])) #AB
                    elif split_line[0] == 'simpopB_tmp':
                        over_approach.append(float(split_line[2])) #BB
            s11_s10_s01_s00_AABBAB.append(over_approach)
        s11_s10_s01_s00 = np.array([[s11_s10_s01_s00_AABBAB[i][j] for i in [0,1,2,3]] for j in [0,2,1]])
        s11_s10_s01_s00 = s11_s10_s01_s00 / np.vstack(np.array([total_bp, total_bp, total_bp])) #Normalised STDs, ie divide by total bp. Alternative would be divide by average.
        #print(s11_s10_s01_s00)
        
        ##NOTE: no power correction but I still need to save each option
        for power in power_grid:
            for powerCorr in powerCorr_grid:
                sim_results_dict_std[','.join([ITx, Ipx, ITa, Ipa, ITb, Ipb, '%.3f' %(power), '%.3f' %(powerCorr)])] = s11_s10_s01_s00

    for filetest_name in list(sim_filebase_set_merged):
        #print(filetest_name)
        #Work out what the parameters are and record the output
        ITx = filetest_name.split('ITx=')[1].split('-')[0]
        ITa = filetest_name.split('ITa=')[1].split('-')[0]
        ITb = filetest_name.split('ITb=')[1].split('-')[0]
        Ipx = filetest_name.split('Ipx=')[1].split('-')[0]
        Ipa = filetest_name.split('Ipa=')[1].split('-')[0]
        Ipb = filetest_name.split('Ipb=')[1].split('-')[0]

        with open(filetest_name) as f_popovermerge:
            s11_s10_s01_s00 = np.array(f_popovermerge.readline().split(','), dtype = float)

        if type(total_bp) == type(None):
            total_bp = np.sum(s11_s10_s01_s00)
        else:
            assert total_bp == np.sum(s11_s10_s01_s00)

        s11_s10_s01_s00 = s11_s10_s01_s00 / float(total_bp)
        
        for power in power_grid:
            for powerCorr in powerCorr_grid:
                sim_results_dict_merged[','.join([ITx, Ipx, ITa, Ipa, ITb, Ipb, '%.3f' %(power), '%.3f' %(powerCorr)])] = modify_11_10_01_00_on_binom(s11_s10_s01_s00, p = 1.0 - power, q = 1.0 - power, rho = powerCorr)


    for filetest_name in list(sim_filebase_set_haps):
        #print(filetest_name)
        #Work out what the parameters are and record the output
        ITx = filetest_name.split('ITx=')[1].split('-')[0]
        ITa = filetest_name.split('ITa=')[1].split('-')[0]
        ITb = filetest_name.split('ITb=')[1].split('-')[0]
        Ipx = filetest_name.split('Ipx=')[1].split('-')[0]
        Ipa = filetest_name.split('Ipa=')[1].split('-')[0]
        Ipb = filetest_name.split('Ipb=')[1].split('-')[0]
        hap_data = []
        with open(filetest_name) as f_hapsdat:
            headers = f_hapsdat.readline()
            for line in f_hapsdat:
                split_line = line[0:-1].split(',')
                hap_data.append(np.array(split_line, dtype = float))
        hap_data = np.array(hap_data)
        pop_split = int(len(hap_data) / 2)

        #Haplotype length is based on a weighted average, so individuals with 0 haplotypes (which are returned as 0 haplotype length too) are skipped
        
        av_haps_A = np.average(hap_data[::,1][0:pop_split])
        av_haps_B = np.average(hap_data[::,1][pop_split:])
        av_hapnum_A = np.average(hap_data[::,2][0:pop_split])
        av_hapnum_B = np.average(hap_data[::,2][pop_split:])
        av_haplen_A = np.sum(hap_data[::,3][0:pop_split] * hap_data[::,2][0:pop_split]) / float(av_hapnum_A * pop_split)
        av_haplen_B = np.sum(hap_data[::,3][pop_split:] * hap_data[::,2][pop_split:]) / float(av_hapnum_B * pop_split)
        
        av_haps_A_internal = np.average(hap_data[::,4][0:pop_split])
        av_haps_B_internal = np.average(hap_data[::,4][pop_split:])
        av_hapnum_A_internal = np.average(hap_data[::,5][0:pop_split])
        av_hapnum_B_internal = np.average(hap_data[::,5][pop_split:])
        av_haplen_A_internal = np.sum(hap_data[::,6][0:pop_split] * hap_data[::,5][0:pop_split]) / float(av_hapnum_A_internal * pop_split)
        av_haplen_B_internal = np.sum(hap_data[::,6][pop_split:] * hap_data[::,5][pop_split:]) / float(av_hapnum_B_internal  * pop_split)

        av_haps_A_internal_over10kb = np.average(hap_data[::,7][0:pop_split])
        av_haps_B_internal_over10kb = np.average(hap_data[::,7][pop_split:])
        av_hapnum_A_internal_over10kb = np.average(hap_data[::,8][0:pop_split])
        av_hapnum_B_internal_over10kb = np.average(hap_data[::,8][pop_split:])
        av_haplen_A_internal_over10kb = np.sum(hap_data[::,9][0:pop_split] * hap_data[::,8][0:pop_split]) / float(av_hapnum_A_internal_over10kb * pop_split) 
        av_haplen_B_internal_over10kb = np.sum(hap_data[::,9][pop_split:] * hap_data[::,8][pop_split:]) / float(av_hapnum_B_internal_over10kb * pop_split) 

        av_haps_A_internal_over20kb = np.average(hap_data[::,10][0:pop_split])
        av_haps_B_internal_over20kb = np.average(hap_data[::,10][pop_split:])
        av_hapnum_A_internal_over20kb = np.average(hap_data[::,11][0:pop_split])
        av_hapnum_B_internal_over20kb = np.average(hap_data[::,11][pop_split:])
        av_haplen_A_internal_over20kb = np.sum(hap_data[::,12][0:pop_split] * hap_data[::,11][0:pop_split]) / float(av_hapnum_A_internal_over20kb * pop_split)
        av_haplen_B_internal_over20kb = np.sum(hap_data[::,12][pop_split:] * hap_data[::,11][pop_split:]) / float(av_hapnum_B_internal_over20kb * pop_split)

        ##NOTE: I do make a power correction for the total and the number only, not the lengths.
        for power in power_grid:
            for powerCorr in powerCorr_grid:
                sim_results_dict_haps[','.join([ITx, Ipx, ITa, Ipa, ITb, Ipb, '%.3f' %(power), '%.3f' %(powerCorr)])] = np.array([av_haps_A * power, av_hapnum_A * power, av_haplen_A,
                                                                                                                                  av_haps_A_internal * power, av_hapnum_A_internal * power, av_haplen_A_internal,
                                                                                                                                  av_haps_A_internal_over10kb * power, av_hapnum_A_internal_over10kb * power, av_haplen_A_internal_over10kb,
                                                                                                                                  av_haps_A_internal_over20kb * power, av_hapnum_A_internal_over20kb * power, av_haplen_A_internal_over20kb,
                                                                                                                                  av_haps_B * power, av_hapnum_B * power, av_haplen_B,
                                                                                                                                  av_haps_B_internal * power, av_hapnum_B_internal * power, av_haplen_B_internal,
                                                                                                                                  av_haps_B_internal_over10kb * power, av_hapnum_B_internal_over10kb * power, av_haplen_B_internal_over10kb,
                                                                                                                                  av_haps_B_internal_over20kb * power, av_hapnum_B_internal_over20kb * power, av_haplen_B_internal_over20kb])
    print("Total states Av: %d Std: %d PopMerge: %d Haps: %d , writing" %(len(sim_results_dict_av), len(sim_results_dict_std), len(sim_results_dict_merged), len(sim_results_dict_haps)))
    
    with open(gridfile_out, 'wb') as f_out:
        headers = ['I_Tx', 'I_prop_x', 'I_Ta',
                   'I_prop_a', 'I_Tb', 'I_prop_b',
                   'power', 'powerCorrRho',
                   's11_AA', 's10_AA', 's01_AA', 's00_AA',
                   's11_BB', 's10_BB', 's01_BB', 's00_BB',
                   's11_AB', 's10_AB', 's01_AB', 's00_AB',
                   's11_AA_std', 's10_AA_std', 's01_AA_std', 's00_AA_std',
                   's11_BB_std', 's10_BB_std', 's01_BB_std', 's00_BB_std',
                   's11_AB_std', 's10_AB_std', 's01_AB_std', 's00_AB_std',
                   's11_AB_popmerge', 's10_AB_popmerge', 's01_AB_popmerge', 's00_AB_popmerge',
                   'av_haps_A', 'av_hapnum_A', 'av_haplen_A',
                   'av_haps_A_internal', 'av_hapnum_A_internal', 'av_haplen_A_internal',
                   'av_haps_A_internal_over10kb', 'av_hapnum_A_internal_over10kb', 'av_haplen_A_internal_over10kb',
                   'av_haps_A_internal_over20kb', 'av_hapnum_A_internal_over20kb', 'av_haplen_A_internal_over20kb',
                   'av_haps_B', 'av_hapnum_B', 'av_haplen_B',
                   'av_haps_B_internal', 'av_hapnum_B_internal', 'av_haplen_B_internal',
                   'av_haps_B_internal_over10kb', 'av_hapnum_B_internal_over10kb', 'av_haplen_B_internal_over10kb',
                   'av_haps_B_internal_over20kb', 'av_hapnum_B_internal_over20kb', 'av_haplen_B_internal_over20kb'
                   ]
        f_out.write((','.join(headers) + '\n').encode())
        for key in np.sort(list(sim_results_dict_av.keys())):
            if (key in sim_results_dict_std.keys()) and (key in sim_results_dict_haps.keys()) and (key in sim_results_dict_merged.keys()):
                line_to_write = ','.join([key, ','.join([','.join(['%.6f' %(j) for j in i ]) for i in sim_results_dict_av[key]])])
                line_to_write = ','.join([line_to_write, ','.join([','.join(['%.6f' %(j) for j in i ]) for i in sim_results_dict_std[key]])])
                line_to_write = ','.join([line_to_write, ','.join([','.join(['%.6f' %(j) for j in sim_results_dict_merged[key]])])])
                line_to_write = ','.join([line_to_write, ','.join([','.join(['%.6f' %(j) for j in sim_results_dict_haps[key]])])])
                #print(line_to_write)
                f_out.write((line_to_write + '\n').encode())
    return len(sim_filebase_set_av)

def assess_dist_grid_obs(f_ingrid, f_outgrid, f_indata_base, groups = ['continent_sasia', 'continent_easia'], metric = 'kl', fitmode = 1):
    return None

def add_derived_parameter_to_fitgrid(f_ingrid, f_outgrid, new_param = 'I_prop_total', target_columns = ['I_prop_x', 'I_prop_a', 'I_prop_b'], operation = lambda x : np.sum(x)):
    with open(f_ingrid, 'rb') as f_in:
        headers = f_in.readline().decode()[0:-1].split(',')
        if type(new_param) == list:
            headers.append(new_param[0])
        else:
            headers.append(new_param)
        target_columns_idx = []
        for targ in target_columns:
            target_columns_idx.append(np.where(np.array(headers) == targ)[0][0])
        lines_to_write = []
        for line in f_in:
            split_line = line.decode()[0:-1].split(',')
            vals = []
            for targ_idx in range(len(target_columns)):
                vals.append(float(split_line[target_columns_idx[targ_idx]]))
            vals = np.array(vals)
            lines_to_write.append(split_line + ['%.6f' %(operation(vals))])
    with open(f_outgrid, 'wb') as f_out:
        f_out.write((','.join(headers) + '\n').encode())
        for line in lines_to_write:
            f_out.write((','.join(line) + '\n').encode())
    return None

def correl_binom(p, q, rho):
    # From https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables
    # Two correlated binomial random numbers to generate the 'false negative' mask
    # p00 has no effect, p01 [01 -> 00; 11 -> 10], p10 [10 -> 00; 11 -> 01], p11 [10 -> 00, 01 -> 00, 11 -> 11]
    a = ((1.0 - p) * (1.0 - q)) + (rho * np.sqrt(p * q * (1.0 - p) * (1.0 - q)))
    p00 = a
    p10 = 1.0 - q - a
    p01 = 1.0 - p - a
    p11 = a + p + q - 1.0
    return p11, p10, p01, p00

def modify_11_10_01_00_on_binom(pobs = np.array([0.25,0.25,0.25,0.25]), p = 0.5, q = 0.5, rho = 1.0):
    # Convert an s11_s10_s01_s00 array into an array with some level of false negatives and some correlation between false negatives
    # p and q are the probability of false negatives in the first and second position respectively. Usually assumed equal.
    pmask = correl_binom(p, q, rho)
    pobs_mod = (np.array([-pobs[0], -pobs[1], -pobs[2], (pobs[1] + pobs[2] + pobs[0])]) * pmask[0]) #m11
    pobs_mod += (np.array([-pobs[0], -pobs[1], pobs[0], pobs[1]]) * pmask[1]) #m10
    pobs_mod += (np.array([-pobs[0], pobs[0], -pobs[2], pobs[2]]) * pmask[2]) #m01
    return pobs + pobs_mod



if args.calc_dist == True:
    raise RuntimeError("Error: distance calculation is not implemented in this version yet")
elif (args.calc_dist == True) and len(args.populations) != 2:
    raise RuntimeError("Error: calculation of distance metric requested but two populations not provided")
elif len(args.derived_parameter_targs) != len(args.derived_parameter_operations) or len(args.derived_parameter_names) != len(args.derived_parameter_operations):
    raise RuntimeError("Error: make sure the number of derived parameter lists equals the number of operations stated and the number of names given.")
elif args.outfile_base[-4:] == '.csv':
    raise RuntimeError("Error: leave off the .csv for the output file, this should be a base file name.")
elif len(args.populations) not in [0,2]:
    raise RuntimeError("Error: the populations should be either 0 for no fitting, or 2 for fitting between two populations.")
elif args.grid_mode not in ['av_only', 'all']:
    raise RuntimeError("Error: the grid_mode should only be av_only (for average haplotype sharing patters only) or all (for all summarised simulated output).")
elif args.archmode not in ['deni', 'nean']:
    raise RuntimeError("Error: the archmode should only be deni (for Denisovan introgression) or nean (for Neanderthal introgression).")
else:
    ## Run the operation
    if (args.archmode == 'deni' and args.arch_power != 0.29) or (args.archmode == 'nean' and args.arch_power != 0.39):
        print(RuntimeWarning("Warning: arch power is not set to the expected value for Jacobs et al 2019 haplotypes. Check?"))

    ## Process the operations into lambda functions
    lambda_list = []
    for operation in args.derived_parameter_operations:
        if 'lambda' in operation[0]:
            lambda_list.append(create_lambda_with_globals(operation[0]))
        elif operation[0] == 'sum':
            lambda_list.append(lambda a : np.sum(a))
        elif operation[0] == 'subtract':
            lambda_list.append(lambda a : a[0] - a[1])
        elif operation[0] == 'multiply':
            lambda_list.append(lambda a : np.product(a))
        elif operation[0] == 'divide':
            lambda_list.append(lambda a : float(a[0]) / float(a[1]))

    ## Run the functions

    if os.path.exists(os.path.dirname(args.outfile_base)) == False:
        try:
            os.mkdir(os.path.dirname(args.outfile_base))
        except FileExistsError:
            pass
    

    current_processing_file = args.outfile_base + '.csv'

    #This creates the csv grid of ABC iterations generated by ABC_sample_IIM_multiAdmix()
    if args.grid_mode == 'av_only':
        num_states = construct_grid_from_ABCdir_avOnly(ABCdir = args.simfolder, gridfile_out = args.outfile_base + '.csv', power_min_max_res = [args.arch_power,args.arch_power,1], powerCorr_min_max_res = args.arch_corr)
    elif args.grid_mode == 'all':
        num_states = construct_grid_from_ABCdir_allData(ABCdir = args.simfolder, gridfile_out = args.outfile_base + '.csv', power_min_max_res = [args.arch_power,args.arch_power,1], powerCorr_min_max_res = args.arch_corr)

    if len(args.derived_parameter_targs) > 0:
        print("Calaculting derived statistics")
        #This adds derived parameters like the total introgression proportion to the csv grid.
        for op_idx in range(len(args.derived_parameter_targs)):
            a = add_derived_parameter_to_fitgrid(f_ingrid = args.outfile_base + '.csv', f_outgrid = args.outfile_base + '.csv', target_columns = args.derived_parameter_targs[op_idx], operation = lambda_list[op_idx], new_param = args.derived_parameter_names[op_idx])

    if args.calc_dist == True:
        #This adds a distance metric to the csv grid
        a = assess_dist_grid_obs(f_ingrid = args.outfile_base + '.csv', f_outgrid = args.outfile_base + '.dist.csv', f_indata_base = args.obsfolder, groups = args.populations, metric = 'kl', fitmode = 1)



