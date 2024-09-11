##Python 3.x

import msprime

import numpy as np
import copy
import pylab as pl
import time
import gzip
import os
import re
import tempfile
import subprocess
import scipy.cluster.hierarchy as hc

import seaborn as sns


SIM_DIR = os.getcwd() + '/raw/'
COMB_DIR = os.getcwd() + '/combined/'


##This conducts the haplotype analysis, simulating introgression overlap based on the Malaspinas et al 2016 model and then comparing results to observations.


###***###   Block Overlap    ###***###

### Creating the simulations to check chunk overlap patterns in Malaspinas-like simulated data ###

### Define a Malaspinas 2016 model, and a function that combines the raw simulation outputs. Functions from msprime_malaspinas2016_out_of_africa_one_wave_py3_v7.

def msprime_malaspinas2016_out_of_africa_one_wave(s_af = 1, s_gh = 1, s_eu = 1, s_ea = 1, s_au = 1, s_deni_i = 1, s_deni_s = 1, s_nean_i = 1, s_nean_s = 1, random_seed = None, num_replicates = 1, length = 1e4):
    """
    This is implemented following Vitor Sousa's fastsimcoal code and powerpoint slide.
    It describes a single OOA (generating the 'ghost' or E African population) followed by a split into Australia and Eurasia, with archaic introgression from Deni and Nean only.
    There are slight differences between the code and the SI; but I cannot really understand these.
    The code is apparently using a generation time of 25 and mutation rate of 1.4e-8, as compared to the SI reported values of 29 generations and 1e-8.
    These are the parameters that Vitor recommends however, so I should follow them.
    """
    # Order of populations is:
    # 0 W Africa
    # 1 Ghost
    # 2 Europe
    # 3 East Asia
    # 4 Australia
    # 5 Deni (introgressing)
    # 6 Deni (sampled)
    # 7 Nean (introgressing)
    # 8 Nean (sampled)

    # Constants:
    mutation_rate = 1.4e-8 # per generation per site, in Vitor's code. NB that 1.25e-8 with generation time 29 was also explored
    generation_time = 25.0 # years, in Vitor's code. NB that 29 years with mutation rate 1.25e-8 was also used.
    recombination_rate = 1.0e-8 # per generation per site. In Malaspinas no recombination was used as the SFS is independent of it. NB that I could optionally sample the recombination rate.

    # Assumed parameters from Prufer et al 2014
    T_divergence_deni_deni_sample = 11998.0 # 13580.0 in SI, gen time effect?
    T_divergence_nean_nean_sample = 3375.0 # 3820.0 in SI
    T_divergence_deni_nean = 15090.0 # 17080.0 in SI
    T_sample_deni_sample = 2058.0 # 2330.0 in SI
    T_sample_nean_sample = 2612.0 # 2957.0 in SI
    inbreeding_nean_sample = 0.125 # Probably the best way to do this is to sample three individuals in t_sample_nean_sample - 2, mate to generate half siblings at t - 1, and then mate to generate half-sib offspring at t.
    
    # Set out the maximum likelihood values of the various parameters given Vitor's code (Nb some differences vs Table S07.3 and S07.5 due ot gentime and mutrate)
    # Ne values are number of diploids. Haploid numbers are given in Vitor's code.
    N_ancestral_archaichuman = 32671 / 2.0 # 18296.0 from SI
    N_deni_altai = 5083 / 2.0 # 2846.0 from SI
    N_unsampled_archaics = 13249 / 2.0 # 7419.0 from SI
    N_nean_altai = 826 / 2.0 # 463.0 from SI
    N_west_africa_yoruba = 48433 / 2.0 # 27122.0 from SI
    N_ghost = 8516 / 2.0 # 4769.0 from SI
    N_europe_sardinia = 6962 / 2.0 # 3899.0 from SI
    N_east_asia_han = 9025 / 2.0 # 5054.0 from SI
    N_aboriginal_australia = 8834 / 2.0 # 4947.0 from SI
    N_ancestral_eurasia = 12971 / 2.0 # 7264.0 from SI
    N_ancestral_humans = 41563 / 2.0 # 23275.0 from SI
    
    N_bottleneck_australia = 243 / 2.0 # 136.0 from SI
    N_bottleneck_eurasia = 2331 / 2.0 # 1305.0 from SI
    N_bottleneck_ooa = 1394 / 2.0 # 781.0 from SI

    #for i in [N_ancestral_archaichuman, N_deni_altai, N_unsampled_archaics, N_nean_altai, N_west_africa_yoruba, N_ghost, N_europe_sardinia, N_east_asia_han, N_aboriginal_australia, N_ancestral_eurasia, N_ancestral_humans, N_bottleneck_australia, N_bottleneck_eurasia, N_bottleneck_ooa]:
    #    #Same order as Vitor's powerpoint
    #    print(i * 2, i)
    
    ## Times are provided in generations in Vitor's code
    T_diverge_human_archaic = 20225.0 #656908.0 / generation_time
    T_diverge_westafrica_ghost = 3916.0 #127192.0 / generation_time
    T_divergence_aboriginal_australians = 1784.0 #57944.0 / generation_time
    T_divergence_eurasia = 1758.0 #57100.0 / generation_time
    T_divergence_european_east_asia = 1293.0 #41997.0 / generation_time
    
    T_bottleneck_ooa = 2119.0 #72041.0 / generation_time
    T_bottleneck_ooa_recovery = 2218.0 #72041.0 / generation_time
    T_bottleneck_eurasia = 1659.0 #This bottleneck stops when the population merges into human ghost
    T_bottleneck_aboriginal_australia = 1685.0 #This bottleneck stops when the population merges into human ghost
    T_admixture_nean_ooa = 1853.0 #60185.0 / generation_time
    T_admixture_nean_eurasia = 1566.0 #50864.0 / generation_time
    T_admixture_nean_east_asia = 883.0 #28680.0 / generation_time
    T_admixture_nean_aboriginal_australia = 1412.0 #45862.0 / generation_time
    T_admixture_deni_aboriginal_australia = 1353.0 #43945.0 / generation_time

    #for i in [T_diverge_westafrica_ghost, T_diverge_human_archaic, T_bottleneck_ooa, T_divergence_aboriginal_australians, T_divergence_eurasia, T_divergence_european_east_asia, T_admixture_nean_eurasia, T_admixture_deni_aboriginal_australia, T_admixture_nean_aboriginal_australia, T_admixture_nean_east_asia, T_admixture_nean_ooa]:
    #    #Same order as Vitor's powerpoint
    #    print(i, i * generation_time)

    ## Introgression proportions
    I_nean_ooa = 0.0235332367
    I_nean_eurasia = 0.0114789187
    I_nean_east_asia = 0.00227784946
    I_nean_aboriginal_australia = 0.00216918517
    I_deni_aboriginal_australia = 0.0397628185
    
    ## Migration rates. In the SI, these were given in units 2Nm; they are symmetrical, with differences reflecting the 2Nm scaling. Rates from powerpoint.

    m_ghost_to_west_africa = 0.000178702636605257
    m_west_africa_to_ghost = 0.000178702636605257
    m_europe_to_ghost = 0.000441607111456954
    m_ghost_to_europe = 0.000441607111456954
    m_ghost_to_ancestral_eurasia = 0.000441607111456954
    m_ancestral_eurasia_to_ghost = 0.000441607111456954
    m_east_asia_to_europe = 3.13921250492658e-05
    m_europe_to_east_asia = 3.13921250492658e-05
    m_australia_to_east_asia = 5.72006968965064e-05
    m_east_asia_to_australia = 5.72006968965064e-05
    m_ancestral_australia_to_ancestral_eurasia = 0.000571694298303235
    m_ancestral_eurasia_to_ancestral_australia = 0.000571694298303235
    
    ## population_configurations_contemporary_sample is a dummy configuration used for the DemographyDebugger
    
    population_configurations_contemporary_sample = [
        msprime.PopulationConfiguration(
            sample_size=s_af, initial_size=N_west_africa_yoruba),
        msprime.PopulationConfiguration(
            sample_size=s_gh, initial_size=N_ghost),
        msprime.PopulationConfiguration(
            sample_size=s_eu, initial_size=N_europe_sardinia),
        msprime.PopulationConfiguration(
            sample_size=s_ea, initial_size=N_east_asia_han),
        msprime.PopulationConfiguration(
            sample_size=s_au, initial_size=N_aboriginal_australia),
        msprime.PopulationConfiguration(
            sample_size=s_deni_i, initial_size=N_unsampled_archaics),
        msprime.PopulationConfiguration(
            sample_size=s_deni_s, initial_size=N_deni_altai),
        msprime.PopulationConfiguration(
            sample_size=s_nean_i, initial_size=N_unsampled_archaics),
        msprime.PopulationConfiguration(
            sample_size=s_nean_s, initial_size=N_nean_altai)
    ]
    
    ## population_configurations is like a population_configurations_contemporary_sample without sample sizes specified to allow for aDNA sampling
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_west_africa_yoruba),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_ghost),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_europe_sardinia),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_east_asia_han),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_aboriginal_australia),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_unsampled_archaics),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_deni_altai),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_unsampled_archaics),
        msprime.PopulationConfiguration(
            sample_size=None, initial_size=N_nean_altai)
    ]

    ## Set up the migration matrix
    
    migration_matrix = [
        [       0.0,            m_ghost_to_west_africa,         0.0,                0.0,                    0.0,                0.0, 0.0, 0.0, 0.0],
        [m_west_africa_to_ghost,        0.0,            m_europe_to_ghost,          0.0,                    0.0,                0.0, 0.0, 0.0, 0.0],
        [       0.0,            m_ghost_to_europe,              0.0,        m_east_asia_to_europe,          0.0,                0.0, 0.0, 0.0, 0.0],
        [       0.0,                    0.0,            m_europe_to_east_asia,      0.0,            m_australia_to_east_asia,   0.0, 0.0, 0.0, 0.0],
        [       0.0,                    0.0,                    0.0,        m_east_asia_to_australia,       0.0,                0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ]

    #High or low migration rates for testing
    #migration_matrix = [[0.01 for i in range(9)] for j in range(9)]
    #for i in range(9):
    #    migration_matrix[i][i] = 0

    demographic_events = [
        # Event 1: Neanderthal introgression into East Asia
        # fastsimcoal2: 883 7 3 0.00227784946 1 0 0 #First event. Neanderthal introgression into East Asia with proportion 0.227% at time 22067
        msprime.MassMigration(
            time = T_admixture_nean_east_asia, source=3, destination=7, proportion=I_nean_east_asia),
        
        # Event 2: Divergence of Europe and East Asia.
        # fastsimcoal2: 1293 7 6 1 1.86299034 0 1 #Second event. Asia merges into Europe at time 32333 years, and population size becomes 12971 (1.86*previous size of Europe). Migration of Europe and ghost doesn't change in terms of prob of lineage movement. Migration of Eurasian and Australia starts.
        # Set migration rates of Pop 3 to 0 and initiate contact between Pop 2 (Ancestral Eurasia) and Pop 4 (Australia)
        # In Vitor's code:
        # 1. We have a population size change to 12971/2.
        # 2. The migration rate between 3 (Asia) and 2 (Europe) of course stops.
        # 3. The migration between 2 (now Eurasia) and 1 (ghost) remains constant, but here needs to be scaled to population size.
        # 4. The migration between 2 (now Eurasia) and 4 (Australia) is set to about 10x the previous migration between E Asia and Australia; and here needs to be scaled to population size too.
        msprime.MassMigration(
            time = T_divergence_european_east_asia, source=3, destination=2, proportion=1.0),
        msprime.PopulationParametersChange(
            time = T_divergence_european_east_asia, initial_size = N_ancestral_eurasia, population_id = 2),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = 0.0, matrix_index = (3, 2)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = 0.0, matrix_index = (2, 3)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = 0.0, matrix_index = (3, 4)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = 0.0, matrix_index = (4, 3)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = m_ancestral_australia_to_ancestral_eurasia, matrix_index = (2, 4)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = m_ancestral_eurasia_to_ancestral_australia, matrix_index = (4, 2)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = m_ghost_to_ancestral_eurasia, matrix_index = (2, 1)),
        msprime.MigrationRateChange(
            time = T_divergence_european_east_asia, rate = m_ancestral_eurasia_to_ghost, matrix_index = (1, 2)),
        # For safety kill off Pop 3 just pastward (nb not allowed to set to 0)
        msprime.PopulationParametersChange(
            time = T_divergence_european_east_asia + 1, initial_size = 0.00001, population_id = 3),
        
        # Event 3: Denisovan introgression into Australia
        # fastsimcoal2: 1353 8 1 0.0397628185 1 0 keep #Third event. Denisovan intrgression into Australia, with proportion 3.9% at time 33817.
        msprime.MassMigration(
            time = T_admixture_deni_aboriginal_australia, source=4, destination=5, proportion=I_deni_aboriginal_australia),
        
        # Event 4: Neanderthal introgression into Australia.
        # fastsimcoal2: 1412 8 3 0.00216918517 1 0 keep #Fourth event, Neanderthal introgression into Australia, with proportion 0.22% at time 39161
        msprime.MassMigration(
            time = T_admixture_nean_aboriginal_australia, source=4, destination=7, proportion=I_nean_aboriginal_australia),
        
        # Event 5: Neanderthal introgression into Eurasia
        # fastsimcoal2: 1566 6 3 0.0114789187 1 0 keep #Fifth event, Neanderthal introgression into Eurasia with proportion 1.15% at time 39161
        msprime.MassMigration(
            time = T_admixture_nean_eurasia, source=2, destination=7, proportion=I_nean_eurasia),
        
        # Event 6: Start of bottleneck in Eurasia
        # fastsimcoal2: 1659 6 6 0 0.179737565 0 2 #Sixth event, bottleneck of Eurasian to 0.179*Ne = 2331 / 2 = 1166
        # This is (back) a 100 generation bottleneck then divergence. Migration stops at this time in Vitor's code.
        msprime.PopulationParametersChange(
            time = T_bottleneck_eurasia, initial_size = N_bottleneck_eurasia, population_id = 2),
        msprime.MigrationRateChange(time = T_bottleneck_eurasia + 0.001, rate = 0.0),

        # Event 7: Start of bottleneck in Australia
        # fastsimcoal2: 1685 8 8 0 0.02755685 0 2 #Seventh event, bottleneck of Australian to 0.02755685 * Ne = 243.4/2 = 122
        # This is (back) a 100 generation bottleneck then divergence.
        msprime.PopulationParametersChange(
            time = T_bottleneck_aboriginal_australia, initial_size = N_bottleneck_australia, population_id = 4),

        # Event 8: Divergence of Eurasia from Ghost.
        # fastsimcoal2: 1758 6 5 1 1 0 2 #Eighth event, Eurasia merges into Ghost and ghost stops migrating (nb after this is was migrating w/ Africa). Population size of ghost doesn't change. Time is 43960 
        msprime.MassMigration(
            time = T_divergence_eurasia, source = 2, destination = 1, proportion = 1.0),
        # Set Eurasia (Pop 2) to 0.0001 for safety
        msprime.PopulationParametersChange(
            time = T_divergence_eurasia + 1, initial_size = 0.00001, population_id = 2),

        # Event 9: Divergence of Australia from Ghost.
        # fastsimcoal2: 1784 8 5 1 1 0 2 #Ninth event, Oceania merges into Ghost and ghost stops migrating (nb after this is was migrating w/ Africa). Population size of ghost doesn't change. Time is 44603
        msprime.MassMigration(
            time = T_divergence_aboriginal_australians, source=4, destination=1, proportion=1.0),
        # Set Australia (Pop 4) to 0.0001 for safety
        msprime.PopulationParametersChange(
            time = T_divergence_aboriginal_australians + 1, initial_size = 0.00001, population_id = 4),


        # Event 10: Neanderthal introgression into OOA
        # fastsimcoal2: 1853 5 3 0.0235332367 1 0 2 #Tenth event, Neanderthal introgression into Ghost, 2.3%. Time is 46335.
        msprime.MassMigration(
            time = T_admixture_nean_ooa, source = 1, destination = 7, proportion = I_nean_ooa),

        # Event 11: OOA bottleneck lasting 100 generations in Ghost (E Africa) population.
        # fastsimcoal2: 2119 5 5 0 0.163722984 0 2 #Eleventh event, Ghost population size changes to 0.163722984 * 8516 = 1394/2 = 697. Time is 52975
        msprime.PopulationParametersChange(
            time = T_bottleneck_ooa, initial_size = N_bottleneck_ooa, population_id = 1),

        # Event 12: OOA bottleneck recovery
        # fastsimcoal2: 2218 5 5 0 6.10787793 0 2 #Twelth event, Ghost population size increases to 697 *  4257 (diploid) = 8516 (haploid).
        msprime.PopulationParametersChange(
            time = T_bottleneck_ooa_recovery, initial_size = N_ghost, population_id = 1),

        # Event 13: Divergence of Nean introgressed and Nean sample
        # fastsimcoal2: 3375 2 3 1 1 0 2 #13th event, Sampled Nean merges with Intro Nean at time 84375. Population size remains at 13249/2 = 6625
        msprime.MassMigration(
            time = T_divergence_nean_nean_sample, source = 8, destination = 7, proportion = 1.0),
        # Set the size of the Neanderthal Ancestors to the unsampled archaic size
        msprime.PopulationParametersChange(
            time = T_divergence_nean_nean_sample, initial_size = N_unsampled_archaics, population_id = 7),
        # Set size of Pop 8 (Introgressing Nean) to 0.00001 for safety
        msprime.PopulationParametersChange(
            time = T_divergence_nean_nean_sample + 0.001, initial_size = 0.00001, population_id = 8),
        
        # Event 14: Divergence of Ghost from West Africa
        # fastsimcoal2: 3916 5 4 1 0.858167331 0 2 #14th event, Ghost merges with Africa at time 97889. Population size of Africa is set to 0.858 * 48433 = 41563/2 = 20778
        msprime.MassMigration(
            time = T_diverge_westafrica_ghost, source=1, destination=0, proportion=1.0),
        # Set the population size of Africa to the ancestral humans size.
        msprime.PopulationParametersChange(
            time = T_diverge_westafrica_ghost, initial_size = N_ancestral_humans, population_id = 0),
        # Set the population size of ghost to 0.00001 for safety
        msprime.PopulationParametersChange(
            time = T_diverge_westafrica_ghost + 0.001, initial_size = 0.00001, population_id = 1),
        
        # Event 15: Divergence of Deni introgressed and Deni sample
        # fastsimcoal2: 11998 0 1 1 1 0 2 #15th event, Deni sampled merges with Deni introgressed at time 299950. Population size stays as it is in Pop 1 i.e. the introgressing Deni, 13249/2 = 6625
        msprime.MassMigration(
            time = T_divergence_deni_deni_sample, source = 6, destination = 5, proportion = 1.0),
        # Set the size of the Denisovan Ancestors to the universal unsampled archaic size
        msprime.PopulationParametersChange(
            time = T_divergence_deni_deni_sample, initial_size = N_unsampled_archaics, population_id = 5),
        # Set size of Pop 6 (Introgressing Deni) to 0.00001 for safety
        msprime.PopulationParametersChange(
            time = T_divergence_deni_deni_sample + 0.001, initial_size = 0.00001, population_id = 6),
        
        # Event 16: Divergence of ancestral Deni and ancestral Nean
        # fastsimcoal2: 15090 1 3 1 1 0 2 #16th event, Deni ancestral and Nean ancesral merge with population size staying at 13249/2 = 6625
        msprime.MassMigration(
            time = T_divergence_deni_nean, source = 7, destination = 5, proportion = 1.0),
        # Set the size of the Nean/Deni Ancestor to the universal unsampled archaic size
        msprime.PopulationParametersChange(
            time = T_divergence_deni_nean, initial_size = N_unsampled_archaics, population_id = 5),
        # Set size of Pop 7 (Ancestral Nean) to 0.00001 for safety
        msprime.PopulationParametersChange(
            time = T_divergence_deni_nean + 0.001, initial_size = 0.00001, population_id = 7),

        # Event 17: Divergence of all hominins
        # fastsimcoal2: 20225 4 3 1 2.46597954 0 2 #Human and archaic merge with population size changing to 2.46597954 * 13249 = 32671 / 2 = 16335
        msprime.MassMigration(
            time = T_diverge_human_archaic, source = 5, destination = 0, proportion = 1.0),
        # Set size of Pop 0 (Nean/Deni/Human Ancestor) to the NAHA size, 16335
        msprime.PopulationParametersChange(
            time = T_diverge_human_archaic, initial_size = N_ancestral_archaichuman, population_id = 0),
        # Set size of Pop 5 (Nean/Deni Ancestor) to 0.00001 for safety
        msprime.PopulationParametersChange(
            time = T_diverge_human_archaic + 0.001, initial_size = 0.00001, population_id = 5),

    ]

    demographic_times = np.array([event.time for event in demographic_events])
    demographic_order = np.argsort(demographic_times)
    demographic_events_sorted = [demographic_events[i] for i in demographic_order]
    
    # Use the demography debugger to print out the demographic history that we have just described.

    N_A = 1.0
    dp = msprime.DemographyDebugger(
        Ne=N_A,
        population_configurations=population_configurations_contemporary_sample,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events_sorted)
    #dp.print_history() # Display the demographic history for debugging

    # NOTE: debugger has to use population_configurations_contemporary_sample as it doens't allow the use of a 'samples' parameter to specify population configuration.
    # The simulation needs to use 'population_configurations' and 'samples' to include aDNA

    # Build the full list of samples to allow for aDNA. Note the order is from pop 0 to pop 8, with sample idx in order
    samples = []
    for pop_idx in range(9):
        # Sampling is assumed to be at time 0, except for the sampled aDNA archaics.
        sample_time = T_sample_deni_sample if pop_idx == 6 else T_sample_nean_sample if pop_idx == 8 else 0.0
        num_samples = [s_af, s_gh, s_eu, s_ea, s_au, s_deni_i, s_deni_s, s_nean_i, s_nean_s][pop_idx]
        samples.extend([msprime.Sample(pop_idx, sample_time) for sample_idx in range(num_samples)])
    #print(samples)
    simulation = msprime.simulate(sample_size = None, #Don't specify, use samples
                                  Ne = N_A,
                                  length = length,
                                  recombination_rate = recombination_rate, # per generation per site
                                  recombination_map = None, # sample this from the HapMap map
                                  mutation_rate = mutation_rate, # per generation per site
                                  population_configurations = population_configurations, #Don't specify, use samples
                                  migration_matrix = migration_matrix,
                                  demographic_events = demographic_events_sorted,
                                  samples = samples,
                                  random_seed = random_seed,
                                  num_replicates = num_replicates,
                                  record_migrations = True)

    #This is the simulation to return.
    
    return simulation



### Functions to work with the migration tables to get introgression sources

def retrieve_chunks_from_simulation_pulse(replicate, target_samples, mrca_samples, target_population = 4, archaic_introgressing_populations = [5, 7], introgressing_time_limits = [[0, 13580], [0, 3820]]):
    """
    This assumes that each chunk either migrated from an archaic (introgressed) or didn't, and that it can only have migrated from one archaic.
    It therefore only works with one-way introgression pulses or events, and doens't allow for e.g. Nean/Deni contact.
    For each target_sample in target_samples, it determines which chunks introgressed from each archaic_sample.
    It also labels the MRCA and location of each chunk with each archaic, and each other target.
    """
    
    introgression_dict = {target:[] for target in target_samples}
    # For each individual, I will have N arrays, for not introgressive, introgression_A, introgression_B... Each array has coalescent information vs targets and archaics:
    # [from, to, coal_arch1_time, coal_arch1_location, ..., coal_targ1_time, coal_targ2_location, ...]
    
    # Go through the migration table.
    # Get every migration that is consistent with introgression, i.e. has a source in the introgression population.
    # That migration may split later on into multiple chunks, i.e. contain multiple coalescent histories
    # So after finding the chunks that are consistent with introgression I need to go through every tree and if:
    # The node corresponding to the chunk is in the tree
    # The relevant interval of that node overlaps with the tree
    # The node has target samples as leaves in the tree
    
    # then it is introgressed and saved as such, with information about it's coalescence also recovered.

    migration_table = replicate.dump_tables().migrations

    #print(migration_table)
    
    migration_nodes = migration_table.node
    migration_sources = migration_table.source
    migration_destinations = migration_table.dest
    migration_times = migration_table.time
    migration_left = migration_table.left
    migration_right = migration_table.right
    start_time = time.time()
    
    for introgression_idx in range(len(archaic_introgressing_populations)):
        introgressed_by_sample = [[] for sample in target_samples]

        introgression_mask = migration_destinations == archaic_introgressing_populations[introgression_idx]

        if introgressing_time_limits[introgression_idx][1] - introgressing_time_limits[introgression_idx][0] > 5:
            print("Introgression bounds are %d and %d. Consider reducing introgressing time limits if possible? With large limits there is a greater potential to incorrectly identify chunks as introgressing. This is especially the case if there are mergers between archaic populations in the specified times, which can lead introgressing chunks to register as coming from multiple archaic populations due to (real) ancestry sharing between them." %(introgressing_time_limits[introgression_idx][0], introgressing_time_limits[introgression_idx][1]))
        within_time_mask = (migration_times[introgression_mask] >= introgressing_time_limits[introgression_idx][0]) * (migration_times[introgression_mask] < introgressing_time_limits[introgression_idx][1])
        
        introgressing_nodes = migration_nodes[introgression_mask][within_time_mask]
        introgressing_left = migration_left[introgression_mask][within_time_mask]
        introgressing_right = migration_right[introgression_mask][within_time_mask]
        
        # A single node can have > 2 children, but only 2 for each chunk within it.
        # I need to ask if the migration being referred to is ancestral to a target sequence.

        #print(target_samples)
        for tree in replicate.trees(tracked_samples = target_samples):
            #print(tree.interval)
            tree_introgressed = set()
            for node_idx in range(len(introgressing_nodes)):
                node = introgressing_nodes[node_idx]
                if (tree.interval[0] >= introgressing_left[node_idx] and tree.interval[0] < introgressing_right[node_idx]) or (tree.interval[1] > introgressing_left[node_idx] and tree.interval[1] <= introgressing_right[node_idx]) or (tree.interval[0] < introgressing_left[node_idx] and tree.interval[1] > introgressing_right[node_idx]) or (tree.interval[0] == introgressing_left[node_idx] and tree.interval[1] == introgressing_right[node_idx]):
                    # If the tree interval is covered, partially or fully, by the migrating node then it is introgressed.
                    if tree.num_tracked_samples(node) > 0:
                        #The introgressing tree has at least one of our target_samples as a leaf; i.e. it is a relevant introgressing tree.
                        tree_leaves_node = set(tree.get_leaves(node))
                        for target in range(len(target_samples)):
                            if (target_samples[target] in tree_leaves_node) and (target not in tree_introgressed):
                                coal_data = list(tree.interval)
                                for mrca_sample in mrca_samples:
                                    if type(mrca_sample) == type(None):
                                        # No MRCA supplied; e.g. a population that would normally be MRCA-ed is not sampled.
                                        # I propose setting these to -1 as they will compress appropriately, unlike np.nan
                                        tmrca = -1
                                        pmrca = -1
                                    elif type(mrca_sample) == int:
                                        mrca_node = tree.get_mrca(target_samples[target], mrca_sample)
                                        tmrca = tree.get_time(mrca_node)
                                        pmrca = tree.get_population(mrca_node)
                                    else:
                                        #mrca_sample is a list or array; I want the *first* common ancestry with any individual in that array
                                        tmrca = np.infty
                                        pmrca = np.nan
                                        for mrca_option in mrca_sample:
                                            mrca_node = tree.get_mrca(target_samples[target], mrca_option)
                                            tmrca_option = tree.get_time(mrca_node)
                                            if tmrca_option < tmrca:
                                                tmrca = copy.deepcopy(tmrca_option)
                                                pmrca = tree.get_population(mrca_node)
                                    coal_data = coal_data + [tmrca, pmrca]
                                introgressed_by_sample[target].append(coal_data)
                                tree_introgressed.add(target)
            
        for target in range(len(target_samples)):
            introgression_dict[target_samples[target]].append(np.array(introgressed_by_sample[target]))
    
    mid_time = time.time()
    # Insert those regions that are not introgressed as index 0
    for target in target_samples:
        target_non_introgressed = []
        archaic_to_check = []
        for introgression_idx in range(len(archaic_introgressing_populations)):
            if len(introgression_dict[target][introgression_idx]) > 0:
                archaic_to_check.append(introgression_idx)
        for tree in replicate.trees():
            non_introgressed = True
            for introgression_idx in archaic_to_check:
                if tree.interval in introgression_dict[target][introgression_idx][::,0:2]:
                    non_introgressed = False
                    break
            if non_introgressed == True:
                coal_data = list(tree.interval)
                for mrca_sample in mrca_samples:
                    if type(mrca_sample) == type(None):
                        # No MRCA supplied; e.g. a population that would normally be MRCA-ed is not sampled.
                        # I propose setting these to -1 as they will compress appropriately, unlike np.nan
                        tmrca = -1
                        pmrca = -1
                    elif type(mrca_sample) == int:
                        mrca_node = tree.get_mrca(target, mrca_sample)
                        tmrca = tree.get_time(mrca_node)
                        pmrca = tree.get_population(mrca_node)
                    else:
                        #mrca_sample is a list or array; I want the *first* common ancestry with any individual in that array
                        tmrca = np.infty
                        pmrca = np.nan
                        for mrca_option in mrca_sample:
                            mrca_node = tree.get_mrca(target, mrca_option)
                            tmrca_option = tree.get_time(mrca_node)
                            if tmrca_option < tmrca:
                                tmrca = copy.deepcopy(tmrca_option)
                                pmrca = tree.get_population(mrca_node)
                    coal_data = coal_data + [tmrca, pmrca]
                target_non_introgressed.append(coal_data)
        introgression_dict[target].insert(0, target_non_introgressed)
        for i in range(len(introgression_dict[target])):
            introgression_dict[target][i] = np.array(introgression_dict[target][i])
    #print(mid_time - start_time, time.time() - mid_time)
    return introgression_dict

def compress_introgression_dict_nocoal(introgression_dict):
    """
    As compress_introgression_dict. In many cases (e.g. visualisation, some data analysis) the coalescent history isn't needed.
    In such cases, I compress but do not keep the coalescent history.
    """
    new_introgression_dict = {}
    for key in introgression_dict.keys():
        new_introgression_dict[key] = []
        for arch in range(len(introgression_dict[key])):
            new_introgression_dict[key].append([])
            curr_idx = 0
            chain_idx = 0
            for val in range(len(introgression_dict[key][arch]) - 1):
                if (introgression_dict[key][arch][val][1] != introgression_dict[key][arch][val + 1][0]):
                    # And start the new one
                    old_chain = [introgression_dict[key][arch][curr_idx][0], introgression_dict[key][arch][chain_idx][1]]
                    new_introgression_dict[key][-1].append(old_chain)
                    curr_idx = val + 1
                    chain_idx = val + 1
                else:
                    # Old chain
                    chain_idx += 1 #val
            #if chain_idx != curr_idx: #GSJ: 25/03/20 Bug here that missed out the (usually small) final chunk, if it was different to the penultimate one
            if len(introgression_dict[key][arch]) > 0:
                final_chain = [introgression_dict[key][arch][curr_idx][0], introgression_dict[key][arch][chain_idx][1]]
                new_introgression_dict[key][-1].append(final_chain)
            new_introgression_dict[key][arch] = np.array(new_introgression_dict[key][arch])
    return new_introgression_dict


def save_mismatch_keep_ind(outfile, mismatch_list):
    #This accepts a series of from, to, mismatch from a single simulation; the mismatch value is optional.
    #It doesn't differentiate which replicate, etc. but does expect ind information to be recorded
    #It saves output in [ind or 0],from,to,mismatch_1, mismatch_2, ... format
    to_write = []
    windows_saved = False
    out_open = gzip.open if outfile[-3:] == '.gz' else open
    for mismatch in mismatch_list:
        # This is [[from_ind0, to_ind0, mis_ind0], ...], [from_ind1, to_ind1, mis_ind1], ...], ...]
        # Convert to [[0, from_ind0, to_ind0, mis_ind0], [0, from_ind0, to_ind0, mis_ind0], [1, from_ind1, to_ind1, mis_ind1], ...]
        mis_by_ind = []
        for ind in range(len(mismatch)):
            for chunk in mismatch[ind]:
                if len(chunk) > 2:
                    #ie there are mismatch values to write
                    mis_by_ind.append([ind, chunk[0], chunk[1], chunk[2]])
                else:
                    mis_by_ind.append([ind, chunk[0], chunk[1], []])
        
        if len(mis_by_ind) != 0:
            mismatch_sort = sorted(mis_by_ind, key = lambda x : (x[1], x[2], x[0]))
            for chunk in range(len(mismatch_sort)):
                if windows_saved == False:
                    to_write.append([mismatch_sort[chunk][0], mismatch_sort[chunk][1], mismatch_sort[chunk][2]])
                else:
                    assert mismatch_sort[chunk][0] == to_write[chunk][0]
                    assert mismatch_sort[chunk][1] == to_write[chunk][1]
                    assert mismatch_sort[chunk][2] == to_write[chunk][2]
                to_write[chunk].append(mismatch_sort[chunk][3])
            windows_saved = True
        else:
            pass
    with out_open(outfile, 'wb') as f:
        for line in to_write:
            try:
                line_str = ','.join(['%d'%(line[0]),'%d'%(line[1]),'%d'%(line[2])] + [('%.3f' %(i) if (np.isnan(i) == True or abs(np.round(i,0) - float(i)) > 1e-8) else '%d' %(np.round(i,0))) for i in np.array(line[3:]).flatten()]) + '\n'
            except:
                print(line)
                raise RuntimeError()
            f.write(line_str.encode())
    return None


### Function to conduct the simulations. NB this version just runs sims and extracts haps.

def experiment_Malaspinas2016_hapoverlap(targ_pop_list = ['au'], individuals_list = [1], num_africans = 1, chromosomes = 5, chrom_len = 1e8, save_chroms = False):
    """
    This is an experiment that applies the Malaspinas et al 2016 model to simulate and study introgression chunks.
    
    """
    assert num_africans >= 1
    assert 'af' not in targ_pop_list
    ref_den_nea_her_mis_den_nea_her_rate = [[[], [], [], [], [], [], [], []] for i in range(4)]
    s_vals = {'af':num_africans, 'gh':0, 'eu':0, 'ea':0, 'au':0}
    for targ_pop_idx in range(len(targ_pop_list)):
        s_vals[targ_pop_list[targ_pop_idx]] = individuals_list[targ_pop_idx]
    target_population = [{'gh':1, 'eu':2, 'ea':3, 'au':4}[targ_pop] for targ_pop in targ_pop_list]
    non_af_individuals = np.sum(individuals_list)
    
    sim = msprime_malaspinas2016_out_of_africa_one_wave(s_af = s_vals['af'],
                                                                         s_gh = s_vals['gh'],
                                                                         s_eu = s_vals['eu'],
                                                                         s_ea = s_vals['ea'],
                                                                         s_au = s_vals['au'],
                                                                         s_deni_i = 0,
                                                                         s_deni_s = 2,
                                                                         s_nean_i = 0,
                                                                         s_nean_s = 2,
                                                                         random_seed = None,
                                                                         num_replicates = chromosomes,
                                                                         length = chrom_len,
                                                                         )
    sim_start_time = time.time()
    # Order is Deni introgression to Australian, Nean introgression to Australian, Her introgression to Australian
    introgressing_time_limits = [[1352.0, 1354.0], [882.0, 1854.0]]
    
    for rep in range(chromosomes):
        rep_start_time = time.time()
        print("Chromosome replicate %d, starting time %f" %(rep, rep_start_time))
        # Reserve the file name for local work
        if save_chroms != False:
            # Reserving file names to save all introgressed chunks. In some versions I save mismatch of these chunks etc.
            for source in ['human', 'deni','nean']:
                f = gzip.open(save_chroms %(rep) + '.blocks.source-%s.BED.gz' %(source), 'wb')
                f.close()

        chrom = next(sim)
        
        """Pre-calculation of introgression chunks and hom-masked geno/pos for replicate"""
        #It is much faster to calculate the introgression_dict and compress it with individuals at once rather than case by case in the mismatch function
        introgression_dict = retrieve_chunks_from_simulation_pulse(chrom,
                                                                   target_samples = range(s_vals['af'], s_vals['af'] + non_af_individuals),
                                                                   mrca_samples = [],
                                                                   target_population = target_population,
                                                                   archaic_introgressing_populations = [5, 7],
                                                                   introgressing_time_limits = introgressing_time_limits)
        data_compress_no_coal = compress_introgression_dict_nocoal(introgression_dict)
        
        for source_idx in range(3):
            source = ['human', 'deni', 'nean'][source_idx]
            if save_chroms != False:
                save_mismatch_keep_ind(outfile = save_chroms %(rep) + '.blocks.source-%s.BED.gz' %(source), mismatch_list = [[[[j[0], j[1]] for j in data_compress_no_coal[i][source_idx]] for i in np.sort(list(data_compress_no_coal.keys()))]])
            
    print("Completed simulation of the Malaspinas et al 2016 model in time %f, total %d inds and %d African outgroup, %.1fbp" %(time.time() - sim_start_time, non_af_individuals, num_africans, chromosomes * chrom_len))
    
    return ref_den_nea_her_mis_den_nea_her_rate



### Function to combine the output of simulations

def combine_mismatch_as_individuals(basefile_mismatch, basefile_out, reps_to_combine, ind_names):
    # This combines all mismatch (or blocks or mrca files), splitting them into individuals.
    # ind_names is a list of VCF ids i.e. ['msp_2', 'msp_3'] if IDxs 0 and 1 in the block files refer to VCF IDs msp_2 and msp_3.
    ind_names_chrom = []
    for ind in ind_names:
        ind_names_chrom.append(ind + '.1')
        ind_names_chrom.append(ind + '.2')
    inds_dict = {ind_name : [] for ind_name in ind_names_chrom}
    open_mismatch = gzip.open if basefile_mismatch[-3:] == '.gz' else open
    for rep in reps_to_combine:
        with open_mismatch(basefile_mismatch %(str(rep)), 'rb') as f_in_base:
            for line in f_in_base:
                split_line = re.split('[,\t]', line[0:-1].decode())
                ind = ind_names[int(np.floor(int(split_line[0]) / 2.0))] #integer division is important
                chrom_copy = (int(split_line[0]) % 2) + 1
                ind_chrom = ind + '.%d' %(chrom_copy)
                split_line[0] = str(rep)
                inds_dict[ind_chrom].append('\t'.join(split_line) + '\n')
    for ind in ind_names_chrom:
        with open_mismatch(basefile_out %('COMBINED.' + ind), 'wb') as f_out_base:
            for line in inds_dict[ind]:
                f_out_base.write(line.encode())
    return None

def combine_simulation_output_blocks(basefile_sims, basefolder_in, basefolder_out, reps_to_combine, ind_names, vcf_reps_are_chroms = False, vcf_reps_are_chunks = False, vcf_reps_are_chunks_insets = None, rename_contig_and_filter_rep_vcfs = False, blocks = False, sources = []):
    # Combines the data of the blocks files
    for source in sources:
        combine_mismatch_as_individuals(basefile_mismatch = basefolder_in + basefile_sims + '.blocks.source-%s.BED.gz' %(source), basefile_out = basefolder_out + basefile_sims + '.blocks.source-%s.BED.gz' %(source), reps_to_combine = reps_to_combine, ind_names = ind_names)
    return None


# Run the simulations. 10 Mb * 100 = 1 Gb  data.
# I do 25 of each population.

a = experiment_Malaspinas2016_hapoverlap(targ_pop_list = ['eu','ea','au'], individuals_list = [50,50,50], num_africans = 1, chromosomes = 100, chrom_len = 1e7, save_chroms = SIM_DIR + 'malaspinasStandardModelEuEaAu_rep%d')

# Merge the outputs

a = combine_simulation_output_blocks(basefile_sims = 'malaspinasStandardModelEuEaAu_rep%s', basefolder_in = SIM_DIR, basefolder_out = COMB_DIR, reps_to_combine = range(100), ind_names = ['eu_%d' %(i) if i < 25 else 'ea_%d' %(i - 25) if i < 50 else 'au_%d' %(i - 50) for i in range(75)], sources = ['deni', 'nean', 'human'])


###*** Analyse the data ***###

BEDTOOLS_DIR = '/home/guy/PhylogenPrograms/bedtools2/bin/'

PIPE = ' | '
OUT = ' > '
POP_DATA = {}

# Functions reading in sample lists and returning where the simulatins for samples

def readin_populationFileEOL(infile):
    with open(infile, 'rb') as f:
        pop = []
        for line in f:
            lins_de = line.decode()
            ind = lins_de[0:-2] if lins_de[-2:] == '\r\n' else lins_de[0:-1] if (lins_de[-1:] == '\r' or lins_de[-1:] == '\n') else lins_de
            pop.append(ind)
    return pop

def input_sims_MalaspinasEuAsAu_deni():
    #Read in the simulated chunks for all individuals (25 European, 25 Asian, 25 Papuan).
    indir = os.getcwd() + '/combined/'
    block_files = []
    ind_names = []
    individuals = ['eu_%d' %(i) if i < 25 else 'ea_%d' %(i - 25) if i < 50 else 'au_%d' %(i - 50) for i in range(75)]
    for ind in individuals:
        for chrom in [1,2]:
            ind_block_files = []
            ind_names.append(ind + '_%d' %(chrom))
            block_files.append(indir + 'malaspinasStandardModelEuEaAu_repCOMBINED.%s.%d.blocks.source-deni.BED.gz' %(ind, chrom))
    return np.array(block_files), np.array(ind_names)

def input_sims_MalaspinasEuAsAu_nean():
    #Read in the simulated chunks for all individuals (25 European, 25 Asian, 25 Papuan).
    indir = os.getcwd() + '/combined/'
    block_files = []
    ind_names = []
    individuals = ['eu_%d' %(i) if i < 25 else 'ea_%d' %(i - 25) if i < 50 else 'au_%d' %(i - 50) for i in range(75)]
    for ind in individuals:
        for chrom in [1,2]:
            ind_block_files = []
            ind_names.append(ind + '_%d' %(chrom))
            block_files.append(indir + 'malaspinasStandardModelEuEaAu_repCOMBINED.%s.%d.blocks.source-nean.BED.gz' %(ind, chrom))
    return np.array(block_files), np.array(ind_names)

def input_skovHaps_nean():
    #Read in the Skov chunks for all individuals, Nean.
    indir = os.getcwd() + '/SkovHaps/'
    block_files = []
    ind_names = []
    block_files = []
    for bed_file in os.listdir(indir):
        if bed_file[-8:] == 'nean.bed':
            block_files.append(indir + bed_file)
            ind_names.append(bed_file.split('.')[0])
    return np.array(block_files), np.array(ind_names)

def input_skovHaps_deni():
    #Read in the Skov chunks for all individuals, Deni.
    indir = os.getcwd() + '/SkovHaps/'
    block_files = []
    ind_names = []
    block_files = []
    for bed_file in os.listdir(indir):
        if bed_file[-8:] == 'deni.bed':
            block_files.append(indir + bed_file)
            ind_names.append(bed_file.split('.')[0])
    return np.array(block_files), np.array(ind_names)

for i in [os.getcwd() + '/popfiles/' + f for f in os.listdir(os.getcwd() + '/popfiles/')]:
    if os.path.isfile(i):
        try:
            list_name = i.split('list_')[1].split('.txt')[0]
            POP_DATA[list_name] = readin_populationFileEOL(i)
        except:
            pass


## Functions for calculating archaic haplotype overlaps between individuals and groups

def popAssess_pairwisePopOverlap_simple(popbedtarg_a = [], popbedtarg_inds = [], popbedcomp_a = [], popbedcomp_inds = [], pairwise_reciprocal = False, distance_approach = 'jaccard', genome_file_local = None, genome_bed_local = None, genome_mask_local = None):
    """
    Assess the overlap between all individuals in one population and all individuals in a second population.
    If reciprocal is true then assess overlap both ways.
    
    This function is used to look at individual overlap; for example, if you want to calculate on all individuals you just pass all individuals in the popbedtarg_a and popbedcomp_a.
    
    """
    ind_overlaps = {}
    for ind_idx in range(len(popbedtarg_a)):
        ind = popbedtarg_inds[ind_idx]
        print('Comparing %s with %s (excluding ideniticals)' %(ind, ' '.join(popbedcomp_inds)))
        for comp_idx in range(len(popbedcomp_a)):
            comp_ind = popbedcomp_inds[comp_idx]
            print(ind, comp_ind)
            if ind != comp_ind and (ind not in ind_overlaps.keys() or (ind in ind_overlaps.keys() and comp_ind not in ind_overlaps[ind])):
                # Compare the two individuals based on coverage or jaccard distance.
                if distance_approach == 'coverage':
                    #Coverage is 11 / [11 + 10] and is asymmetrical
                    overlap_command = BEDTOOLS_DIR + 'bedtools coverage -a %s -b %s -g %s' %(popbedtarg_a[ind_idx], popbedcomp_a[comp_idx], genome_file_local) + PIPE + "awk '{SUMcover+=$5;SUMtot+=$6}END{if (SUMtot > 0) {print SUMcover/SUMtot} else {print 0}}'"
                    output_transform = lambda a : np.array([float(a)])
                elif distance_approach == 'jaccard':
                    #Jaccard is 11 / [11 + 10 + 01] and so is symmetrical
                    overlap_command = BEDTOOLS_DIR + 'bedtools jaccard -a %s -b %s -g %s' %(popbedtarg_a[ind_idx], popbedcomp_a[comp_idx], genome_file_local)
                    output_transform = lambda a : np.array([float(a.split(b'\n')[1].split(b'\t')[2])])
                elif distance_approach == 'all':
                    
                    overlap_command_11 = BEDTOOLS_DIR + 'bedtools coverage -a %s -b %s -g %s' %(popbedtarg_a[ind_idx], popbedcomp_a[comp_idx], genome_file_local) + PIPE + "awk '{SUMcover+=$5}END{print SUMcover}'"
                    overlap_command_10 = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -g %s' %(popbedtarg_a[ind_idx], popbedcomp_a[comp_idx], genome_file_local) + PIPE + "awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"
                    overlap_command_01 = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s -g %s' %(popbedcomp_a[comp_idx], popbedtarg_a[ind_idx], genome_file_local) + PIPE + "awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"
                    
                    overlap_command_11_10_01_00 = BEDTOOLS_DIR + 'bedtools subtract -a %s -b %s' %(genome_bed_local, genome_mask_local) + PIPE + "awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"
                    result_11_10_01_00 = []
                    bedtools_process_11 = subprocess.Popen(args = overlap_command_11, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                    stdout, stderr = bedtools_process_11.communicate(input = None)
                    result_11_10_01_00.append(0 if stdout == b'\n' else int(stdout))
                    
                    
                    bedtools_process_10 = subprocess.Popen(args = overlap_command_10, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                    stdout, stderr = bedtools_process_10.communicate(input = None)
                    result_11_10_01_00.append(0 if stdout == b'\n' else int(stdout))
                    
                    bedtools_process_01 = subprocess.Popen(args = overlap_command_01, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                    stdout, stderr = bedtools_process_01.communicate(input = None)
                    result_11_10_01_00.append(0 if stdout == b'\n' else int(stdout))
                    
                    bedtools_process_11_10_01_00 = subprocess.Popen(args = overlap_command_11_10_01_00, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                    stdout, stderr = bedtools_process_11_10_01_00.communicate(input = None)
                    result_11_10_01_00.append(int(stdout) - np.sum(result_11_10_01_00))
                    result_11_10_01_00 = np.array(result_11_10_01_00)
                    
                    #print(ind, comp_ind, result_11_10_01_00)
                    
                    try:
                        ind_overlaps[ind][comp_ind] = result_11_10_01_00
                    except KeyError:
                        ind_overlaps[ind] = {}
                        ind_overlaps[ind][comp_ind] = result_11_10_01_00
                if distance_approach != 'all':
                    #print(overlap_command)
                    bedtools_process_coverage = subprocess.Popen(args = overlap_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                    stdout, stderr = bedtools_process_coverage.communicate(input = None)
                    if stderr is not b'':
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
                        overlap_command = BEDTOOLS_DIR + 'bedtools coverage -a %s -b %s -g %s' %(popbedcomp_a[comp_idx], popbedtarg_a[ind_idx], genome_file_local) + PIPE + "awk '{SUMcover+=$5;SUMtot+=$6}END{if (SUMtot > 0) {print SUMcover/SUMtot} else {print 0}}'"
                    elif distance_approach == 'jaccard':
                        # The actual values of the reciprocal should be idenitical? They appear to be very very close but not necessarily quite identical...
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

                        #print(comp_ind, ind, ind_overlaps[comp_ind][ind])
                    if distance_approach != 'all':
                        # print(overlap_command)
                        bedtools_process_coverage = subprocess.Popen(args = overlap_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                        stdout, stderr = bedtools_process_coverage.communicate(input = None)
                        if stderr is not b'':
                            print(overlap_command)
                            print(stdout)
                            raise RuntimeError(stderr)
                        try:
                            ind_overlaps[comp_ind][ind] = output_transform(stdout)
                        except KeyError:
                            ind_overlaps[comp_ind] = {}
                            ind_overlaps[comp_ind][ind] = output_transform(stdout)
    return ind_overlaps

## This function calculates archaic hap overlaps between populations

def multipleAssess_pairwisePopOverlap_simple(func_filelist, outfile, haploid_calls = True, pops = ['sims_all75'], file_IDx_A = 0, distance_approach = 'jaccard', genome_file_local = None, genome_bed_local = None, genome_mask_local = None):
    """
    Ask for overlap between individuals or populations.
    
    In this case, I also want to average over _1 and _2 of individuals.
    """
    assess_files, ind_names = func_filelist()
    
    order = np.sort(np.array(pops))
    all_results = {}
    av_results = {}
    for pop_targ_idx in np.arange(len(order)):
        #go through each
        pop_targ = order[pop_targ_idx]
        print("Calculating overlaps for target population %s" %(pop_targ))
        pop_results = {}
        poptarget_assess = []
        poptarget_inds = []
        for ind in range(len(ind_names)):
            if haploid_calls == True:
                if ind_names[ind][0:-2] in POP_DATA[pop_targ]:
                    poptarget_assess.append(assess_files[ind])
                    poptarget_inds.append(ind_names[ind])
            else:
                if ind_names[ind] in POP_DATA[pop_targ]:
                    poptarget_assess.append(assess_files[ind])
                    poptarget_inds.append(ind_names[ind])
        poptarget_assess = np.array(poptarget_assess)
        poptarget_inds = np.array(poptarget_inds)
        for pop_comp_idx in np.arange(len(order) - pop_targ_idx) + pop_targ_idx:
            pop_comp = order[pop_comp_idx]
            print("Comparing with comparison population %s" %(pop_comp))
            popcomp_assess = []
            popcomp_inds = []
            for ind in range(len(ind_names)):
                if haploid_calls == True:
                    if ind_names[ind][0:-2] in POP_DATA[pop_comp]:
                        popcomp_assess.append(assess_files[ind])
                        popcomp_inds.append(ind_names[ind])
                else:
                    if ind_names[ind] in POP_DATA[pop_comp]:
                        popcomp_assess.append(assess_files[ind])
                        popcomp_inds.append(ind_names[ind])
            popcomp_assess = np.array(popcomp_assess)
            popcomp_inds = np.array(popcomp_inds)
            print(poptarget_assess)
            pop_results_targ = popAssess_pairwisePopOverlap_simple(popbedtarg_a = poptarget_assess,
                                                    popbedtarg_inds = poptarget_inds,
                                                    popbedcomp_a = popcomp_assess,
                                                    popbedcomp_inds = popcomp_inds,
                                                    distance_approach = distance_approach,
                                                    pairwise_reciprocal = True,
                                                    genome_file_local = genome_file_local,
                                                    genome_bed_local = genome_bed_local,
                                                    genome_mask_local = genome_mask_local)
            for ind in pop_results_targ.keys():
                if haploid_calls == True:
                    if ind[0:-2] not in all_results.keys():
                        all_results[ind[0:-2]] = {}
                    for comp in pop_results_targ[ind].keys():
                        if comp[0:-2] not in all_results[ind[0:-2]].keys():
                            all_results[ind[0:-2]][comp[0:-2]] = [pop_results_targ[ind][comp]]
                        else:
                            all_results[ind[0:-2]][comp[0:-2]].append(pop_results_targ[ind][comp])
                else:
                    if ind not in all_results.keys():
                        all_results[ind] = {}
                    for comp in pop_results_targ[ind].keys():
                        if comp not in all_results[ind].keys():
                            all_results[ind][comp] = [pop_results_targ[ind][comp]]
                        else:
                            all_results[ind][comp].append(pop_results_targ[ind][comp])
    print("COMPLETED! Writing.")
    f_list = [outfile] if distance_approach != 'all' else [outfile + '.11.csv', outfile + '.10.csv', outfile + '.01.csv', outfile + '.00.csv']
    for f_out_idx in range(len(f_list)):
        f_out = open(f_list[f_out_idx], 'wb')
        av_results = {}
        for ind in all_results.keys():
            av_results[ind] = {}
            for comp in all_results[ind].keys():
                av_results[ind][comp] = np.nanmean(np.array(all_results[ind][comp])[::,f_out_idx])
        ind_order = np.sort(list(all_results.keys()))
        f_out.write((','.join(['IND'] + list(ind_order)) + '\n').encode())
        for ind in ind_order:
            line_to_write = [ind]
            for comp_ind in ind_order:
                try:
                    line_to_write.append('%.4f' %(av_results[ind][comp_ind]))
                except:
                    line_to_write.append('nan')
            f_out.write((','.join(line_to_write) + '\n').encode())
        f_out.close()
    return all_results


##Analysis step 1:

# Calculate the overlap of simulations using bedtools_archaic_v1.py

a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_sims_MalaspinasEuAsAu_deni, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni', pops = ['sims_all75'], distance_approach = 'all', file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes', genome_bed_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.mask')

a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_sims_MalaspinasEuAsAu_deni, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_coverage_simEuEaAu_simdeni.csv', pops = ['sims_all75'], distance_approach = 'coverage', file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes', genome_bed_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.mask')

a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_sims_MalaspinasEuAsAu_deni, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_jaccard_simEuEaAu_simdeni.csv', pops = ['sims_all75'], distance_approach = 'jaccard', file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes', genome_bed_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.mask')


a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_sims_MalaspinasEuAsAu_nean, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean', pops = ['sims_all75'], distance_approach = 'all', file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes', genome_bed_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.mask')

a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_sims_MalaspinasEuAsAu_nean, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_coverage_simEuEaAu_simnean.csv', pops = ['sims_all75'], distance_approach = 'coverage', file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes', genome_bed_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.mask')

a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_sims_MalaspinasEuAsAu_nean, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_jaccard_simEuEaAu_simnean.csv', pops = ['sims_all75'], distance_approach = 'jaccard', file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes', genome_bed_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/sim_100chunksnum_10Mb.mask')



### Jaccard for Skov too, so I can do a quick comparison with the sims
a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_skovHaps_nean, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_jaccard_all_skovnean.csv', pops = ['subset_skov_all'], distance_approach = 'jaccard', haploid_calls = False, file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/GRCh37_UCSC_autoonly_hg19.chrom.sizes', genome_bed_local = os.getcwd() + '/maskfiles/GRCh37_UCSC_autoonly_hg19.chrom.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/hg19_full_GAP_mask.bed')
a = multipleAssess_pairwisePopOverlap_simple(func_filelist = input_skovHaps_deni, outfile = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_jaccard_all_skovdeni.csv', pops = ['subset_skov_all'], distance_approach = 'jaccard', haploid_calls = False, file_IDx_A = 0, genome_file_local = os.getcwd() + '/maskfiles/GRCh37_UCSC_autoonly_hg19.chrom.sizes', genome_bed_local = os.getcwd() + '/maskfiles/GRCh37_UCSC_autoonly_hg19.chrom.sizes.bed', genome_mask_local = os.getcwd() + '/maskfiles/hg19_full_GAP_mask.bed')


# Average the matrices to populations and get out the StDev. This could be useful to know.

##Analysis step 2:

### Function to average pairwise overlaps into a pop X pop matrix

def average_pairwise_ind_matrix_to_pop(in_matrix_ind = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.11.csv', out_matrix_pop = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.11.pops.%s.csv', pops = ['sims_eu25', 'sims_ea25']):
    """
    Read in a matix of individual pairwise comparisons and write out the mean and the standard deviation of
    pairwise comparisons between a given set of populations.
    """
    matrix_ind = {}
    with open(in_matrix_ind, 'rb') as f_in:
        headers = f_in.readline().decode()[0:-1].split(',')
        for line in f_in:
            split_line = line[0:-1].decode().split(',')
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
                print(pop_targ, pop_comp)
                inds_targ = POP_DATA[pop_targ]
                inds_comp = POP_DATA[pop_comp]
                to_compare = []
                for ind_targ in inds_targ:
                    for ind_comp in inds_comp:
                        to_compare.append(matrix_ind[ind_targ][ind_comp])
                        try:
                            matrix_pop_av[pop_targ][pop_comp] = np.nanmean(to_compare)
                        except KeyError:
                            matrix_pop_av[pop_targ] = {}
                            matrix_pop_av[pop_targ][pop_comp] = np.nanmean(to_compare)
                        try:
                            matrix_pop_std[pop_targ][pop_comp] = np.nanstd(to_compare)
                        except KeyError:
                            matrix_pop_std[pop_targ] = {}
                            matrix_pop_std[pop_targ][pop_comp] = np.nanstd(to_compare)
        with open(out_matrix_pop %('av'), 'wb') as f_out:
            headers = ["POP"] + pops
            f_out.write((','.join(headers) + '\n').encode())
            for pop_targ in pops:
                line_to_write = [pop_targ]
                for pop_comp in pops:
                    line_to_write.append('%.4f' %(matrix_pop_av[pop_targ][pop_comp]))
                f_out.write((','.join(line_to_write) + '\n').encode())
        with open(out_matrix_pop %('std'), 'wb') as f_out:
            headers = ["POP"] + pops
            f_out.write((','.join(headers) + '\n').encode())
            for pop_targ in pops:
                line_to_write = [pop_targ]
                for pop_comp in pops:
                    line_to_write.append('%.4f' %(matrix_pop_std[pop_targ][pop_comp]))
                f_out.write((','.join(line_to_write) + '\n').encode())
    return None

def convert_pairwise_to_prop(in_file_base, out_file_base, file_sub_list):
    read_data_all = []
    headers = ''
    for file_sub in file_sub_list:
        read_data = []
        with open(in_file_base %(file_sub), 'rb') as f_in:
            if headers != '':
                headers_new = f_in.readline().decode()
                #print(headers, headers_new)
                assert headers == headers_new
                headers = headers_new
            else:
                headers = f_in.readline().decode()
            for line in f_in:
                read_data.append(np.array(line.decode().strip().split(',')[1:], dtype = float))
        read_data_all.append(read_data)
    headers = headers.strip().split(',')
    num_pops = len(read_data)
    for file_sub_idx in range(len(file_sub_list)):
        file_sub = file_sub_list[file_sub_idx]
        with open(out_file_base %(file_sub), 'wb') as f_out:
            f_out.write((','.join(headers) + '\n').encode())
            for i in range(num_pops):
                line_to_write = [headers[i + 1]]
                for j in range(num_pops):
                    counts = [read_data_all[k][i][j] for k in range(len(file_sub_list))]
                    prop = counts[file_sub_idx] / np.sum(counts)
                    line_to_write.append('%.6f' %(prop))
                f_out.write((','.join(line_to_write) + '\n').encode())
    return None


for archaic in ['deni', 'nean']:
    for target in ['pairwiseIndOver_coverage_allNonSSAf_%sHCSS35unique.csv', 'pairwiseIndOver_jaccard_allNonSSAf_%sHCSS35unique.csv', 'pairwiseIndOver_all_allNonSSAf_%sHCSS35unique.01.csv', 'pairwiseIndOver_all_allNonSSAf_%sHCSS35unique.10.csv', 'pairwiseIndOver_all_allNonSSAf_%sHCSS35unique.11.csv', 'pairwiseIndOver_all_allNonSSAf_%sHCSS35unique.00.csv']:
        dat_dir = os.getcwd() + '/pairwiseIndOver/'
        dat_file = dat_dir + target %(archaic)
        dat_out = dat_file.replace('.csv', '.pops.%s.csv')
        average_pairwise_ind_matrix_to_pop(in_matrix_ind = dat_file, out_matrix_pop = dat_out, pops = ['continent_easia', 'continent_europe', 'continent_sasia', 'subset_papuaExclFrancoisUVBaining'])
    convert_pairwise_to_prop(in_file_base = os.getcwd() + '/pairwiseIndOver/' + 'pairwiseIndOver_all_allNonSSAf_%sHCSS35unique.%s.pops.av.csv' %(archaic, '%s'), out_file_base = os.getcwd() + '/pairwiseIndOver/' + 'pairwiseIndOver_all_allNonSSAf_%sHCSS35unique.%s.pops.prop.av.csv' %(archaic, '%s'), file_sub_list = ['11', '10', '01', '00'])



for archaic in ['nean', 'deni']:
    for target in ['pairwiseIndOver_coverage_simEuEaAu_sim%s.csv', 'pairwiseIndOver_jaccard_simEuEaAu_sim%s.csv', 'pairwiseIndOver_all_simEuEaAu_sim%s.01.csv', 'pairwiseIndOver_all_simEuEaAu_sim%s.10.csv', 'pairwiseIndOver_all_simEuEaAu_sim%s.11.csv', 'pairwiseIndOver_all_simEuEaAu_sim%s.00.csv']:
        dat_dir = os.getcwd() + '/pairwiseIndOver/'
        dat_file = dat_dir + target %(archaic)
        dat_out = dat_file.replace('.csv', '.pops.%s.csv')
        average_pairwise_ind_matrix_to_pop(in_matrix_ind = dat_file, out_matrix_pop = dat_out, pops = ['sims_eu25', 'sims_ea25', 'sims_au25'])
    convert_pairwise_to_prop(in_file_base = os.getcwd() + '/pairwiseIndOver/' + 'pairwiseIndOver_all_simEuEaAu_sim%s.%s.pops.av.csv' %(archaic, '%s'), out_file_base = os.getcwd() + '/pairwiseIndOver/' + 'pairwiseIndOver_all_simEuEaAu_sim%s.%s.pops.prop.av.csv' %(archaic, '%s'), file_sub_list = ['11', '10', '01', '00'])


### And for Skov
for archaic in ['deni', 'nean']:
    for target in ['pairwiseIndOver_jaccard_all_skov%s.csv']:
        dat_dir = os.getcwd() + '/pairwiseIndOver/'
        dat_file = dat_dir + target %(archaic)
        dat_out = dat_file.replace('.csv', '.pops.%s.csv')
        if archaic == 'nean':
            average_pairwise_ind_matrix_to_pop(in_matrix_ind = dat_file, out_matrix_pop = dat_out, pops = ['subset_skov_easia', 'subset_skov_sasia', 'subset_skov_weurasia'])
        else:
            average_pairwise_ind_matrix_to_pop(in_matrix_ind = dat_file, out_matrix_pop = dat_out, pops = ['subset_skov_easia', 'subset_skov_sasia'])


##Analysis step 3:

#Exploring the impact of power on these results.

### Estimated power to detect Deni is 0.29 and estimated power to detect Nean is 0.39.
### Just consider power correlation of 0, 0.5, 1.0 to get a broad sense of things.

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
    #print(pobs, np.sum(pobs))
    #print(pmask)
    #print(p,q,rho)
    #print(pobs_mod)
    #print(pobs + pobs_mod)
    #print('\n')
    return pobs + pobs_mod


def recalc_jaccard_from_average_11_10_01_00_files(file_base, file_out_base, p = 0.5, q = 0.5, rho = 1.0):
    #file_base is the input base of the 11/10/01/00 matrices, with one %s substitution for the '11' etc strings
    #file_out_base is the ouput base with three %f substitutions for the p, q and rho strings
    site_count_all = []
    for i in ['11', '10', '01', '00']:
        site_counts = []
        with open(file_base %(i), 'rb') as f_in:
            header_line = f_in.readline()
            pop_header = header_line.decode().strip().split(',')[1:]
            for line in f_in:
                av_props = np.array(line.decode().strip().split(',')[1:])
                site_counts.append(av_props)
        site_count_all.append(site_counts)
    #site_count_all is [[11[(popA,popA), (popA,popB), ...],[(popB,popA), (popB,popB), ...], ...], [10[(popA,popA), ...], ...], etc.]
    site_count_all = np.array(site_count_all)
    jaccard_matrix = np.zeros((len(pop_header), len(pop_header)))
    for i in range(len(pop_header)):
        for j in range(len(pop_header)):
            site_counts_pair_mod = modify_11_10_01_00_on_binom(pobs = np.array([site_count_all[0][i][j], site_count_all[1][i][j], site_count_all[2][i][j], site_count_all[3][i][j]], dtype = float), p = p, q = q, rho = rho)
            jaccard_pair_mod = site_counts_pair_mod[0] / np.sum([site_counts_pair_mod[1], site_counts_pair_mod[2]])
            jaccard_matrix[i][j] = jaccard_pair_mod
    with open(file_out_base %(p, q, rho), 'wb') as f_out:
        print(pop_header)
        header = ['POP'] + pop_header
        f_out.write((','.join(header) + '\n').encode())
        for line_idx in range(len(jaccard_matrix)):
            line_to_write = ','.join(['%.4f' %(i) for i in jaccard_matrix[line_idx]]) + '\n'
            f_out.write(line_to_write.encode())
    return None

a = recalc_jaccard_from_average_11_10_01_00_files(file_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.%s.pops.prop.av.csv', file_out_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.jaccard.pow_p%.2f_q%.2f_rho%.2f.pops.prop.av.csv', p = 0.29, q = 0.29, rho = 0.0)
a = recalc_jaccard_from_average_11_10_01_00_files(file_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.%s.pops.prop.av.csv', file_out_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.jaccard.pow_p%.2f_q%.2f_rho%.2f.pops.prop.av.csv', p = 0.29, q = 0.29, rho = 0.5)
a = recalc_jaccard_from_average_11_10_01_00_files(file_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.%s.pops.prop.av.csv', file_out_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simdeni.jaccard.pow_p%.2f_q%.2f_rho%.2f.pops.prop.av.csv', p = 0.29, q = 0.29, rho = 1.0)

a = recalc_jaccard_from_average_11_10_01_00_files(file_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean.%s.pops.prop.av.csv', file_out_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean.jaccard.pow_p%.2f_q%.2f_rho%.2f.pops.prop.av.csv', p = 0.39, q = 0.39, rho = 0.0)
a = recalc_jaccard_from_average_11_10_01_00_files(file_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean.%s.pops.prop.av.csv', file_out_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean.jaccard.pow_p%.2f_q%.2f_rho%.2f.pops.prop.av.csv', p = 0.39, q = 0.39, rho = 0.5)
a = recalc_jaccard_from_average_11_10_01_00_files(file_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean.%s.pops.prop.av.csv', file_out_base = os.getcwd() + '/pairwiseIndOver/pairwiseIndOver_all_simEuEaAu_simnean.jaccard.pow_p%.2f_q%.2f_rho%.2f.pops.prop.av.csv', p = 0.39, q = 0.39, rho = 1.0)


##Analysis step 4: Plotting

### Visualisations of cluster plots ###

def figure_pairwiseIndOverlap_cluster(pairwisePopOverlapFile, pairwisePopOverlapFile_comparisonList = [], pairwisePopOverlapFile_comparisonFunction = lambda x, y : x - y, symmetry_method = 'average', data_mode = 'populations', mask_pops = False, colour_limits = [None, None], cmap = 'rocket'):
    #This data shows to what extent chunks overlap.
    #There is a question as the X and Y trees are only the same if the data is symmetrical. I have a few options of how to do that.
    #I can draw and colour various options of populations, continents, ...
    #For comparisons between Deni and Nean, it might be useful to have e.g. colour_limits = [0.0, 0.5]. Because the max overlap is different for Deni and Nean.
    #The first step is to read in the individuals.

    #pairwisePopOverlapFile_comparisonList and pairwisePopOverlapFile_comparisonFunction let me e.g. calculate the residuals from another function.
    #E.g. lambda x, y : x - y is residuals, lambda x, y, z : x / np.sum([x,y,z]) using 11, 10, 01 give proportion of total pairwise detected introgression...
    
    data = []
    data_inds = []
    data_pops = []

    #Simulated populations
    pops_simEuEaAu = ['sims_eu25', 'sims_ea25', 'sims_au25']
    pop_simEuEaAu_names = ['a) Sim:Europe', 'b) Sim:East Asia', 'c) Sim:Australia']
    #Europe/E Asia/Papua
    pop_EuEaPap_names = ['a) West Eurasia', 'b) East Asia', 'c) Papua']
    pops_EuEaPap_25 = ['subset_europe_rnd25', 'subset_easia_rnd25', 'subset_papuaNoUVNG_rnd25']
    #Continents, don't overlap
    continents_exclFrancoisUVBaining = ['subset_papuaExclFrancoisUVBaining', 'continent_e_isea', 'continent_w_isea', 'continent_seasia', 'continent_easia', 'continent_america', 'continent_oceania', 'continent_siberia', 'continent_sasia', 'continent_europe']
    cont_names = ['Papua', 'East ISEA', 'West ISEA', 'SE Asia', 'E Asia', 'America', 'Oceania', 'Siberia', 'S Asia', 'Europe']
    
    if data_mode == 'sims75':
        pops_plot = pops_simEuEaAu
        pops_plot_names = pop_simEuEaAu_names
    elif data_mode == 'eueapap_focus25':
        pops_plot = pops_EuEaPap_25
        pops_plot_names = pop_EuEaPap_names
    elif data_mode == 'continents_exclFrancoisUVBaining':
        pops_plot = continents_exclFrancoisUVBaining
        pops_plot_names = cont_names
    
    
    pops_plot_names = np.array(pops_plot_names)
    pops_plot = np.array(pops_plot)
    
    f_in = open(pairwisePopOverlapFile, 'rb')
    headers = f_in.readline().decode()[0:-1].split(',')
    
    for line in f_in:
        split_line = line.decode()[0:-1].split(',')
        ind = split_line[0]
        ind_pop = 'None'
        for popcheck in pops_plot:
            if ind in POP_DATA[popcheck]:
                if ind_pop == 'None':
                    ind_pop = pops_plot_names[np.where(pops_plot == popcheck)[0][0]]
                else:
                    raise RuntimeError('Ind %s member of multiple populations %s %s' %(ind, ind_pop, popcheck))
            else:
                pass
        data.append([max(float(i), 1e-5) for i in split_line[1:]])
        data_inds.append(ind)
        data_pops.append(ind_pop)
    if mask_pops == True:
        data_pops = np.array(data_pops)
        pop_mask = data_pops != 'None'
        data = np.array(data)[pop_mask]
        data = data[::,pop_mask]
        data_inds = np.array(data_inds)[pop_mask]
        data_pops = data_pops[pop_mask]
    else:
        data_pops = np.array(data_pops)
        data_inds = np.array(data_inds)
        data = np.array(data)

    #colour_set = sns.color_palette('Set1', n_colors = len(np.unique(data_pops[data_pops != 'None'])))
    colour_set = sns.color_palette('Set3', n_colors = len(np.unique(data_pops[data_pops != 'None'])))
    colours = {'None':(0.0,0.0,0.0)}
    for pop in np.sort(np.unique(data_pops)):
        if pop != 'None':
            colours[pop] = colour_set[len(colours)-1]
    rowcol_colours = [colours[pop] for pop in data_pops]
    if symmetry_method == 'average':
        symdat = (data + np.fliplr(np.rot90(data,3))) / 2.0 #Reciprocal average
    elif symmetry_method == 'multiply':
        symdat = (data * np.fliplr(np.rot90(data,3)))
    elif symmetry_method == 'max':
        print(data)
        symdat = np.maximum(data, np.fliplr(np.rot90(data,3)))
        print(symdat)
    elif symmetry_method == 'min':
        print(data)
        symdat = np.minimum(data, np.fliplr(np.rot90(data,3)))
        print(symdat)
    elif symmetry_method == 'againstPops':
        #In this case, I want to summarise each individual as the average distance between each member of
        comp_pop_names = []
        comp_pop_idx = []
        for pop in np.unique(data_pops):
            if pop != 'None':
                comp_pop_idx.append(data_pops == pop)
                comp_pop_names.append(pop)
        new_data = []
        for ind in data:
            for pop in range(len(comp_pop_names)):
                pass
    elif symmetry_method == 'none':
        #Get the average anyway for ordering
        symdat = (data + np.fliplr(np.rot90(data,3))) / 2.0 #Reciprocal average
    linkage_matrix = hc.linkage(y = symdat, method = 'average', optimal_ordering = True)
    
    if symmetry_method != 'none':
        cluster = sns.clustermap(data = symdat, col_linkage = linkage_matrix, row_linkage = linkage_matrix, col_colors = rowcol_colours, xticklabels = data_inds, yticklabels = data_inds, method = 'average', vmin = colour_limits[0], vmax = colour_limits[1], cmap = cmap)
    else:
        cluster = sns.clustermap(data = data, col_linkage = linkage_matrix, row_linkage = linkage_matrix, col_colors = rowcol_colours, xticklabels = data_inds, yticklabels = data_inds, method = 'average', vmin = colour_limits[0], vmax = colour_limits[1], cmap = cmap)
    for x in cluster.ax_heatmap.get_xticklabels():
        x.set_fontsize(3)
        x.set_rotation(90)
    for y in cluster.ax_heatmap.get_yticklabels():
        y.set_fontsize(3)
        y.set_rotation(0)
    leftax = cluster.ax_row_dendrogram.axes
    for item in leftax.collections:
        item.set_visible(False)
    for pop in np.unique(data_pops):
        leftax.bar(0,0,color=colours[pop], label=pop, linewidth=0)
    legend = leftax.legend()
    legend.get_frame().set_visible(False)
    return None


#Simulations:

for archaic in ['nean', 'deni']:
    for target in ['pairwiseIndOver_jaccard_simEuEaAu_sim%s.csv']:
        for approach in ['none']:
            indir = os.getcwd() + '/pairwiseIndOver/'
            min_max = [-0.01, 0.2] if 'jaccard' in target else [-0.01,0.3] if 'coverage' in target else [None, None]
            cmap = pl.cm.get_cmap('bone_r', 1024)
            print(archaic, target, approach)
            a = figure_pairwiseIndOverlap_cluster(pairwisePopOverlapFile = indir + target %(archaic), symmetry_method = approach, data_mode = 'sims75', mask_pops = False, colour_limits = min_max, cmap = cmap)
            outdir = os.getcwd() + '/figs/'
            filename_cluster = outdir + target %(archaic)
            filename_cluster = filename_cluster[0:-4] + '_ClusterSymmetry%s.pdf' %(approach.capitalize())
            pl.savefig(filename_cluster)
            pl.close()

# 25 individual actual data plots

for archaic in ['deni', 'nean']:
    for target in ['pairwiseIndOver_jaccard_allNonSSAf_%sHCSS35unique.csv']:
        for approach in ['none']:
            indir = os.getcwd() + '/pairwiseIndOver/'
            min_max = [-0.01, 0.2] if 'jaccard' in target else [-0.01,0.3] if 'coverage' in target else [None, None]
            print(archaic, target, approach)
            cmap = pl.cm.get_cmap('bone_r', 1024)
            a = figure_pairwiseIndOverlap_cluster(pairwisePopOverlapFile = indir + target %(archaic), symmetry_method = approach, data_mode = 'eueapap_focus25' if archaic == 'nean' else 'eueapap_focus25', mask_pops = True, colour_limits = min_max, cmap = cmap)
            
            outdir = os.getcwd() + '/figs/'
            filename_cluster = outdir + target %(archaic)
            filename_cluster = filename_cluster[0:-4] + '_rnd25' + '_ClusterSymmetry%s.pdf' %(approach.capitalize())
            pl.savefig(filename_cluster)
            pl.clf()
            pl.close()
            
# All individual actual data plots

#NOTE - to get the correct colours add these lines into 
#colours["E Asia"] = sns.color_palette()[0]
#colours["Papua"] = sns.color_palette()[2]
#colours["Europe"] = sns.color_palette()[3]


for archaic in ['deni', 'nean']:
    for target in ['pairwiseIndOver_jaccard_allNonSSAf_%sHCSS35unique.csv']:
        for approach in ['none']:
            indir = os.getcwd() + '/pairwiseIndOver/'
            min_max = [-0.01, 0.2] if 'jaccard' in target else [-0.01,0.3] if 'coverage' in target else [None, None]
            print(archaic, target, approach)
            cmap = pl.cm.get_cmap('bone_r', 1024)
            a = figure_pairwiseIndOverlap_cluster(pairwisePopOverlapFile = indir + target %(archaic), symmetry_method = approach, data_mode = 'continents_exclFrancoisUVBaining', mask_pops = True, colour_limits = min_max, cmap = cmap)
            
            outdir = os.getcwd() + '/figs/'
            filename_cluster = outdir + target %(archaic)
            filename_cluster = filename_cluster[0:-4] + '_allIndsExclFrancoisUVBaining' + '_ClusterSymmetry%s.pdf' %(approach.capitalize())
            pl.savefig(filename_cluster)
            pl.clf()
            pl.close()

##Analysis step 4: Residuals

#GSJ: Added Jan 2023

def figure_pairwiseIndOverlap_cluster_comparison(pairwisePopOverlapFile_simulated75, pairwisePopOverlapFile_sampled75, symmetry_method = 'average', comparison_function = lambda x, y : x / y, colour_limits = [None, None], cmap = 'rocket'):
    #This data shows to what extent chunks overlap.
    #The funciont builds is a trimmed down version of figure_pairwiseIndOverlap_cluster specifically for these residuals
  
    data_sim = []
    data_sim_inds = []
    data_sim_pops = []
    data_samp = []
    data_samp_inds = []
    data_samp_pops = []

    #Simulated populations
    pops_plot_sim = np.array(['sims_eu25', 'sims_ea25', 'sims_au25'])
    pops_plot_sim_names = np.array(['a) Sim:Europe', 'b) Sim:East Asia', 'c) Sim:Australia'])
    #Europe/E Asia/Papua
    pops_plot_EuEaPap = np.array(['subset_europe_rnd25', 'subset_easia_rnd25', 'subset_papuaNoUVNG_rnd25'])
    pops_plot_EuEaPap_names = np.array(['a) West Eurasia', 'b) East Asia', 'c) Papua'])
    
    #Read in simulated data
    f_in_sim = open(pairwisePopOverlapFile_simulated75, 'rb')
    headers_sim = f_in_sim.readline().decode()[0:-1].split(',')
    
    for line in f_in_sim:
        split_line = line.decode()[0:-1].split(',')
        ind = split_line[0]
        ind_pop = 'None'
        for popcheck in pops_plot_sim:
            if ind in POP_DATA[popcheck]:
                if ind_pop == 'None':
                    ind_pop = pops_plot_sim_names[np.where(pops_plot_sim == popcheck)[0][0]]
                else:
                    raise RuntimeError('Ind %s member of multiple populations %s %s' %(ind, ind_pop, popcheck))
            else:
                pass
        data_sim.append([max(float(i), 1e-5) for i in split_line[1:]])
        data_sim_inds.append(ind)
        data_sim_pops.append(ind_pop)

    #Read in the sampled data
    f_in_sample = open(pairwisePopOverlapFile_sampled75, 'rb')
    headers_sample = f_in_sample.readline().decode()[0:-1].split(',')
    
    for line in f_in_sample:
        split_line = line.decode()[0:-1].split(',')
        ind = split_line[0]
        ind_pop = 'None'
        for popcheck in pops_plot_EuEaPap:
            if ind in POP_DATA[popcheck]:
                if ind_pop == 'None':
                    ind_pop = pops_plot_EuEaPap_names[np.where(pops_plot_EuEaPap == popcheck)[0][0]]
                else:
                    raise RuntimeError('Ind %s member of multiple populations %s %s' %(ind, ind_pop, popcheck))
            else:
                pass
        data_samp.append([max(float(i), 1e-5) for i in split_line[1:]])
        data_samp_inds.append(ind)
        data_samp_pops.append(ind_pop)
    
    data_sim_pops = np.array(data_sim_pops)
    data_sim_inds = np.array(data_sim_inds)
    data_sim = np.array(data_sim)

    data_samp_pops = np.array(data_samp_pops)
    data_samp_inds = np.array(data_samp_inds)
    data_samp = np.array(data_samp)

    pop_mask = data_samp_pops != 'None'
    data_samp = np.array(data_samp)[pop_mask]
    data_samp = data_samp[::,pop_mask]
    data_samp_inds = np.array(data_samp_inds)[pop_mask]
    data_samp_pops = data_samp_pops[pop_mask]

    new_order = np.argsort(data_samp_pops)[::-1]
    data_samp = data_samp[new_order]
    data_samp = data_samp[::, new_order]
    data_samp_inds = data_samp_inds[new_order]
    data_samp_pops = data_samp_pops[new_order]
    

    colour_set = sns.color_palette('Set1', n_colors = len(np.unique(data_samp_pops[data_samp_pops != 'None'])))
    #colour_set = sns.color_palette('Set3', n_colors = len(np.unique(data_samp_pops[data_samp_pops != 'None'])))
    colours = {'None':(0.0,0.0,0.0)}
    for pop in np.sort(np.unique(data_samp_pops)):
        if pop != 'None':
            colours[pop] = colour_set[len(colours)-1]
    colours["E Asia"] = sns.color_palette()[0]
    colours["Papua"] = sns.color_palette()[2]
    colours["Europe"] = sns.color_palette()[3]

    rowcol_colours = [colours[pop] for pop in data_samp_pops]
    
    if symmetry_method == 'average':
        symdat_sim = (data_sim + np.fliplr(np.rot90(data_sim,3))) / 2.0 #Reciprocal average
        symdat_samp = (data_samp + np.fliplr(np.rot90(data_samp,3))) / 2.0 #Reciprocal average
        comparison_data = comparison_function(symdat_samp, symdat_sim)
    
    elif symmetry_method == 'none':
        #Get the average anyway for ordering
        symdat_sim = (data_sim + np.fliplr(np.rot90(data_sim,3))) / 2.0 #Reciprocal average
        symdat_samp = (data_samp + np.fliplr(np.rot90(data_samp,3))) / 2.0 #Reciprocal average
        comparison_data = comparison_function(data_samp, data_sim)

    #print(data_samp_pops)
    #print(data_sim_pops)

    #print(data_samp)
    #print(data_sim)
    
    linkage_matrix = hc.linkage(y = symdat_samp, method = 'average', optimal_ordering = True)
    
    cluster = sns.clustermap(data = comparison_data, col_linkage = linkage_matrix, row_linkage = linkage_matrix, col_colors = rowcol_colours, xticklabels = data_samp_inds, yticklabels = data_samp_inds, method = 'average', vmin = colour_limits[0], vmax = colour_limits[1], cmap = cmap)
    

    for x in cluster.ax_heatmap.get_xticklabels():
        x.set_fontsize(3)
        x.set_rotation(90)
    for y in cluster.ax_heatmap.get_yticklabels():
        y.set_fontsize(3)
        y.set_rotation(0)
    leftax = cluster.ax_row_dendrogram.axes
    for item in leftax.collections:
        item.set_visible(False)
    for pop in np.unique(data_samp_pops):
        leftax.bar(0,0,color=colours[pop], label=pop, linewidth=0)
    legend = leftax.legend()
    legend.get_frame().set_visible(False)

for archaic in ['nean', 'deni']:
    for targets in [['pairwiseIndOver_jaccard_simEuEaAu_sim%s.csv', 'pairwiseIndOver_jaccard_allNonSSAf_%sHCSS35unique.csv']]:
        for approach in ['none']:
            indir = os.getcwd() + '/pairwiseIndOver/'
            #min_max = [-0.01, 0.2] if 'jaccard' in targets[0] else [-0.01,0.3] if 'coverage' in targets[0] else [None, None]
            min_max = [-1.0,1.0]
            #min_max = [-0.01,0.1]
            #cmap = pl.cm.get_cmap('bone_r', 1024)
            cmap = pl.cm.get_cmap('BrBG', 1024)
            print(archaic, targets, approach)
            a = figure_pairwiseIndOverlap_cluster_comparison(pairwisePopOverlapFile_simulated75 = indir + targets[0] %(archaic), pairwisePopOverlapFile_sampled75 = indir + targets[1] %(archaic), symmetry_method = approach, comparison_function = lambda x, y : np.log2(x / y), colour_limits = min_max, cmap = cmap)
            outdir = os.getcwd() + '/figs/'
            filename_cluster = outdir + target %(archaic)
            filename_cluster = filename_cluster[0:-4] + '_ClusterSymmetry%s.pdf' %(approach.capitalize())
            pl.savefig(filename_cluster)
            pl.clf()
            pl.close()
