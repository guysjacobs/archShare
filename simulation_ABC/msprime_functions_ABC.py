##Functions to run msprime simulations and process output

import msprime
import numpy as np
import copy
import time
import gzip
import re
import tempfile
import subprocess

PIPE = '|'

def msprime_IIM_model(T0_gens, T1_gens, a, b, c1, c2, N, m1, m2, I_Tx, I_Ta, I_Tb, I_prop_x, I_prop_a, I_prop_b, s_A = 1, s_B = 1, random_seed = None, num_replicates = 1, length = 1e4):
    """
    This is a simple Isolation with Inital Migration model with a pulse of introgression.

    The population:
    1. Has an ancestral size of Na
    2. Splits at T_0 into two populations of size N and Nb, with migration rates m1 and m2 between them
    3. Migration stops at T_1 and the population sizes become Nc1 and Nc2
    4. Sampling occurs at time 0

    There is the option of up to three introgression events from one archaic source.
    These occur with proportion I_prop_x, I_prop_a and I_prop_b and times I_Tx, I_Ta, I_Tb

    I model the archaic source as splitting at 20225 gens, and a Deni effective population size as 7419; both from Malaspinas

    #Example demography of Papua vs East Asia from /IIM/fitting/resupper2/Summary_Table_IIM2only.xlsx
    #T0_gens = 3502, T1_gens = 468, a = 7.37, b = 1.01, c1 = 47.29, c2 = 55.80, N = 2000.0, m1 = 0.00067, m2 = 0.00035

    """
    # Order of populations is:
    # 0 A
    # 1 B
    # 2 Archaic_into_A
    # 3 Archaic_into_B
    # 4 Archaic_into_X

    if I_Tx <= T0_gens:
        raise RuntimeError("I_Tx %.2f (shared introgression) must be more than T0_gens %.2f (split time)" %(I_Tx, T0_gens))
    if I_Tx >= 20225:
        raise RuntimeError("I_Tx %.2f (shared introgression) must be less than 20225 gens (archaic/human split time)" %(I_Tx))
    if I_Ta >= T0_gens:
        raise RuntimeError("I_Ta %.2f (pop A introgression) must be less than T0_gens %.2f (split time)" %(I_Ta, T0_gens))
    if I_Tb >= T0_gens:
        raise RuntimeError("I_Tb %.2f (pop B introgression) must be less than T0_gens %.2f (split time)" %(I_Tb, T0_gens))
    
    # Constants:
    mutation_rate = 1.4e-8 # per generation per site, in Vitor's code. NB that 1.25e-8 with generation time 29 was also explored
    generation_time = 25.0 # years, in Vitor's code. NB that 29 years with mutation rate 1.25e-8 was also used.
    recombination_rate = 1.0e-8 # per generation per site. In Malaspinas no recombination was used as the SFS is independent of it. NB that I could optionally sample the recombination rate.

    T_divergence_arch_hum = 20225.0 # Assumed split at 20225 gens from Malaspinas. Split time of populations 2/3/4 too.
    T_divergence_A_B = float(T0_gens) # Time at which populations A and B diverge
    T_popsizechange_A_B = float(T1_gens) # Time at which populations A and B change size and change migration rate
    
    # Set out the maximum likelihood values of the various parameters given Vitor's code (Nb some differences vs Table S07.3 and S07.5 due ot gentime and mutrate)
    # Ne values are number of diploids. Haploid numbers are given in Vitor's code.
    N_ancestral = N * a # Ancestral size
    N_archaic = 13249 / 2.0 # 7419.0, from Malaspinas
    
    N_A_T0_T1 = N # Population size of A between times T0 and T1
    N_B_T0_T1 = N * b # Population size of B between times T0 and T1
    N_A_T1_0 = N * c1 # Population size of A between times T1 and 0
    N_B_T1_0 = N * c2 # Population size of B between times T1 and 0

    ## Introgression parameters
    T_admixture_arch_X = I_Tx # Generations at which introgression occurs into ancestral population X
    T_admixture_arch_A = I_Ta # Generations at which introgression occurs into population A
    T_admixture_arch_B = I_Tb # Generations at which introgression occurs into population B

    I_admixture_arch_X = I_prop_x # Amount of introgression occuring into ancestral population X
    I_admixture_arch_A = I_prop_a # Amount of introgression occuring into population A
    I_admixture_arch_B = I_prop_b # Amount of introgression occuring into population B

    
    
    ## Migration rates. Backwards in time.
    m_A_to_B = m1
    m_B_to_A = m2
    
    ## population_configurations_contemporary_sample is a dummy configuration used for the DemographyDebugger
    
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=s_A, initial_size=N_A_T1_0),
        msprime.PopulationConfiguration(
            sample_size=s_B, initial_size=N_B_T1_0),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_archaic),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_archaic),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_archaic)
    ]
    
    ## Set up the migration matrix
    
    migration_matrix = [[       0.0,            0.0,            0.0,        0.0,        0.0 ],
                        [       0.0,            0.0,            0.0,        0.0,        0.0 ],
                        [       0.0,            0.0,            0.0,        0.0,        0.0 ],
                        [       0.0,            0.0,            0.0,        0.0,        0.0 ],
                        [       0.0,            0.0,            0.0,        0.0,        0.0 ]]


    demographic_events = [
        # Unsorted. Deal with the human first, then the three possible introgression events.
        # Event 1: Population sizes change and migration starts.

        msprime.PopulationParametersChange(
            time = T_popsizechange_A_B, initial_size = N_A_T0_T1, population_id = 0),
        msprime.PopulationParametersChange(
            time = T_popsizechange_A_B, initial_size = N_B_T0_T1, population_id = 1),
        msprime.MigrationRateChange(
            time = T_popsizechange_A_B + 1, rate = m_A_to_B, matrix_index = (1, 0)),
        msprime.MigrationRateChange(
            time = T_popsizechange_A_B + 1, rate = m_B_to_A, matrix_index = (0, 1)),

        # Event 2: Populations 1 and 2 merger. Migration stops. Population size changes.

        msprime.MassMigration(
            time = T_divergence_A_B, source=1, destination=0, proportion=1.0),
        msprime.PopulationParametersChange(
            time = T_divergence_A_B, initial_size = N_ancestral, population_id = 0),
        msprime.MigrationRateChange(
            time = T_divergence_A_B - 1, rate = 0.0, matrix_index = (0, 1)),
        msprime.MigrationRateChange(
            time = T_divergence_A_B - 1, rate = 0.0, matrix_index = (1, 0)),
        msprime.PopulationParametersChange(
            time = T_divergence_A_B + 1, initial_size = 0.00001, population_id = 1),

        # Event 3: Archaic population 4 merges into archaic population 3
        msprime.MassMigration(
            time = T_divergence_arch_hum - 2, source=4, destination=3, proportion=1.0),
        msprime.PopulationParametersChange(
            time = T_divergence_arch_hum + 1, initial_size = 0.00001, population_id = 4),

        # Event 3: Archaic population 3 merges into archaic population 2
        msprime.MassMigration(
            time = T_divergence_arch_hum - 1, source=3, destination=2, proportion=1.0),
        msprime.PopulationParametersChange(
            time = T_divergence_arch_hum + 1, initial_size = 0.00001, population_id = 3),

        # Event 3: Archaic population merges into the ancestral human population
        msprime.MassMigration(
            time = T_divergence_arch_hum, source=2, destination=0, proportion=1.0),
        msprime.PopulationParametersChange(
            time = T_divergence_arch_hum + 1, initial_size = 0.00001, population_id = 2),

        # Event Intro 1: Introgression into A
        msprime.MassMigration(
            time = T_admixture_arch_A, source=0, destination=2, proportion=I_admixture_arch_A),

        # Event Intro 1: Introgression into B
        msprime.MassMigration(
            time = T_admixture_arch_B, source=1, destination=3, proportion=I_admixture_arch_B),

        # Event Intro 1: Introgression into X
        msprime.MassMigration(
            time = T_admixture_arch_X, source=0, destination=4, proportion=I_admixture_arch_X)

    ]

    demographic_times = np.array([event.time for event in demographic_events])
    demographic_order = np.argsort(demographic_times)
    demographic_events_sorted = [demographic_events[i] for i in demographic_order]
    
    # Use the demography debugger to print out the demographic history that we have just described.

    N_A = 1.0
    dp = msprime.DemographyDebugger(
        Ne=N_A,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events_sorted)
    #dp.print_history() # Display the demographic history for debugging

    # NOTE: debugger has to use population_configurations_contemporary_sample as it doens't allow the use of a 'samples' parameter to specify population configuration.
    # The simulation needs to use 'population_configurations' and 'samples' to include aDNA

    simulation = msprime.simulate(sample_size = None, #Don't specify, use samples s_A + s_B + s_arch
                                  Ne = N_A,
                                  length = length,
                                  recombination_rate = recombination_rate, # per generation per site
                                  recombination_map = None, # sample this from the HapMap map
                                  mutation_rate = mutation_rate, # per generation per site
                                  population_configurations = population_configurations,
                                  migration_matrix = migration_matrix,
                                  demographic_events = demographic_events_sorted,
                                  random_seed = random_seed,
                                  num_replicates = num_replicates,
                                  record_migrations = True)

    #This is the simulation to return.
    
    return simulation


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



### Function to conduct the simulations.


def experiment_IIM_multipleintro(T0_gens = 3502, T1_gens = 468, a = 7.37, b = 1.01, c1 = 47.29, c2 = 55.80, N = 2000.0, m1 = 0.00067, m2 = 0.00035, I_Tx = 3550, I_Ta = 1500, I_Tb = 1500, I_prop_x = 0.03, I_prop_a = 0.02, I_prop_b = 0.01, num_popA = 1, num_popB = 1, chromosomes = 5, chrom_len = 1e8, save_chroms = False):
    """
    This is an experiment involving an IIM simulation with potentially multiple introgression events.
    There are two populations that split, with migration for a period and then isolation.
    There are potential introgression events in the ancestral populatinon and in the two daughter populations.
    I am interested in reporting the BED files of locations of introgression only.

    We do not calculate mismatch of chunks ro return VCFs for this method.

    """
    start_time = time.time()
    s_vals = {'a':num_popA, 'b':num_popB}
    
    sim = msprime_IIM_model(T0_gens = T0_gens,
                            T1_gens = T1_gens,
                            a = b,
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
                            s_A = s_vals['a'],
                            s_B = s_vals['b'],
                            random_seed = None,
                            num_replicates = chromosomes,
                            length = chrom_len)

    # Order is archaic introgression to A, archaic introgression to B, archaic introgression to ancestral
    introgressing_time_limits = [[I_Ta - 1, I_Ta + 1], [I_Tb - 1, I_Tb + 1], [I_Tx - 1, I_Tx + 1]]
    print(introgressing_time_limits)
    
    for rep in range(chromosomes):
        rep_start_time = time.time()
        print("Chromosome replicate %d, starting time %f" %(rep, rep_start_time))
        # Reserve the file name for local work
        sources = ['human','archA', 'archB', 'archX']
        if save_chroms != False:
            # This option saves which blocks are from which introgression source. Note that I try to keep track of which source.
            for source in sources:
                f = gzip.open(save_chroms %(rep) + '.popA.blocks.source-%s.BED.gz' %(source), 'wb')
                f.close()
                f = gzip.open(save_chroms %(rep) + '.popB.blocks.source-%s.BED.gz' %(source), 'wb')
                f.close()
        #print("Presim", time.time())
        chrom = next(sim)
        #print("Postsim", time.time())

        """Calculation of admixture chunks for replicate"""
        #It is much faster to calculate the admixture_dict and compress it with all individuals at once
        time_chunkStart = time.time()
        
        admixture_dict_pop = retrieve_chunks_from_simulation_pulse(chrom,
                                                               target_samples = range(s_vals['a'] + s_vals['b']),
                                                               mrca_samples = [],
                                                               target_population = 1, #Nb target_population doesn't impact this function
                                                               archaic_introgressing_populations = [2,3,4],
                                                               introgressing_time_limits = introgressing_time_limits)
        data_compress_no_coal = compress_introgression_dict_nocoal(admixture_dict_pop)
        data_compress_no_coalA = {x:data_compress_no_coal[x] for x in range(s_vals['a'])}
        data_compress_no_coalB = {x:data_compress_no_coal[x] for x in range(s_vals['a'], s_vals['a'] + s_vals['b'])}
        #print("chunk ops took", time.time() - time_chunkStart)
        
        for source_idx in [0,1,2,3]:
            source = ['human', 'archA', 'archB', 'archX'][source_idx]
            
            if save_chroms != False:
                # This option saves the mismatch of all introgressed chunks, and the mismatch and locations of chunks separately.
                # Note that I still save all individuals in one file. Fine to separate them later I think.
                save_mismatch_keep_ind(outfile = save_chroms %(rep) + '.popA.blocks.source-%s.BED.gz' %(source), mismatch_list = [[[[j[0], j[1]] for j in data_compress_no_coalA[i][source_idx]] for i in np.sort(list(data_compress_no_coalA.keys()))]])
                save_mismatch_keep_ind(outfile = save_chroms %(rep) + '.popB.blocks.source-%s.BED.gz' %(source), mismatch_list = [[[[j[0], j[1]] for j in data_compress_no_coalB[i][source_idx]] for i in np.sort(list(data_compress_no_coalB.keys()))]])
    end_time = time.time()
    print("Completed simulation of IIM model in time %.2f, with parameters: \n[T0,T1] = [%d,%d], \n[a,b,c1,c2,N] = [%.2f,%.2f,%.2f,%.2f,%.2f], \n[m1,m2] = [%.6f,%.6f], \n[I_Tx,I_prop_x] = [%d, %.6f], \n[I_Ta,I_prop_a] = [%d, %.6f], \n[I_Tb,I_prop_b] = [%d, %.6f]" %(end_time - start_time, T0_gens, T1_gens, a, b, c1, c2, N, m1, m2, I_Tx, I_prop_x, I_Ta, I_prop_a, I_Tb, I_prop_b))
    
    return None


###---OUTPUT PROCESSING---###


def merge_beds(outfile, bed_list = [], bedtools_dir = '/home/guy/PhylogenPrograms/bedtools2/bin/', tmp_folder = '/home/guy/tmp/'):
    PIPE = ' | '
    
    files_to_merge = ' '.join(bed_list)
    command_cat_and_merge_list = 'zcat %s' %(files_to_merge) + PIPE + 'sort -k1,1V -k2,2n' + PIPE + bedtools_dir + 'bedtools merge -i -'
    f_cat_and_merge_list = tempfile.NamedTemporaryFile(mode = 'w+b', dir = tmp_folder)
    bedtools_process_cat_and_merge_list = subprocess.Popen(args = command_cat_and_merge_list, stdout = f_cat_and_merge_list, stderr = f_cat_and_merge_list, shell = True)
    bedtools_process_cat_and_merge_list.communicate(input=None)
    f_cat_and_merge_list.flush()
    f_cat_and_merge_list.seek(0)
        
    #And write out
    with open(outfile, 'wb') as f_out_merged:
        for line in f_cat_and_merge_list:
            f_out_merged.write(line)
    f_cat_and_merge_list.close()
    return None
    
def merge_beds_uncompressed(outfile, bed_list = [], bedtools_dir = '/home/guy/PhylogenPrograms/bedtools2/bin/', tmp_folder = '/home/guy/tmp/'):
    PIPE = ' | '
    
    files_to_merge = ' '.join(bed_list)
    command_cat_and_merge_list = 'cat %s' %(files_to_merge) + PIPE + 'sort -k1,1V -k2,2n' + PIPE + bedtools_dir + 'bedtools merge -i -'
    f_cat_and_merge_list = tempfile.NamedTemporaryFile(mode = 'w+b', dir = tmp_folder)
    bedtools_process_cat_and_merge_list = subprocess.Popen(args = command_cat_and_merge_list, stdout = f_cat_and_merge_list, stderr = f_cat_and_merge_list, shell = True)
    bedtools_process_cat_and_merge_list.communicate(input=None)
    f_cat_and_merge_list.flush()
    f_cat_and_merge_list.seek(0)
        
    #And write out
    with open(outfile, 'wb') as f_out_merged:
        for line in f_cat_and_merge_list:
            f_out_merged.write(line)
    f_cat_and_merge_list.close()
    return None

def combine_VCF(basefile_vcf, basefile_out, reps_to_combine, reps_are_chroms_or_chunks = 'chroms'):
    #This combined the VCFs from multiple simulated chromosomes.
    #First, get the headers
    assert basefile_vcf.count('%s') == 1
    open_vcf = gzip.open if basefile_vcf[-3:] == '.gz' else open
    headers = []
    with open_vcf(basefile_vcf %(str(reps_to_combine[0])), 'rb') as f_vcf:
        for line in f_vcf:
            if line[0:1] == '#':
                headers.append(line)
            else:
                break
    contigs = []
    for rep in reps_to_combine:
        try:
            with open_vcf(basefile_vcf %(str(rep)), 'rb') as f_vcf:
                for line in f_vcf:
                    if line[0:2] == '##':
                        if '##contig=' in line:
                            contigs.append(line.replace('<ID=1,', '<ID=%s,' %(str(rep))))
                    elif line[0] == '#':
                        assert line == headers[-1]
                    else:
                        break
        except:
            raise RuntimeError("Couldn't read in rep %s" %(str(rep)))
    with open_vcf(basefile_out %('COMBINED'), 'wb') as f_out_vcf:
        # Write the headers
        for header_id in range(len(headers)):
            header = headers[header_id]
            if '##contig=' in header and reps_are_chroms_or_chunks == 'chroms':
                # Don't write the contig line from the header, which is always 1-idx
                # Unless the reps are chunks of a single chromosome, in which case, keep that single chromosome!
                for contig in contigs:
                    f_out_vcf.write(contig)
            else:
                f_out_vcf.write(header)
        bp_total = 0
        # Write the data
        for rep_idx in range(len(reps_to_combine)):
            rep = reps_to_combine[rep_idx]
            contig = contigs[rep_idx]
            if rep_idx != 0:
                bp_total += int(contigs[rep_idx - 1].split(',length=')[1].split('>')[0])
            with open_vcf(basefile_vcf %(str(rep)), 'rb') as f_vcf_in:
                for line in f_vcf_in:
                    if line[0] == '#':
                        pass
                    else:
                        split_line = line.split('\t')
                        split_line[0] = str(rep) if reps_are_chroms_or_chunks == 'chroms' else split_line[0] if reps_are_chroms_or_chunks == 'chunks' else None
                        split_line[1] = split_line[1] if reps_are_chroms_or_chunks == 'chroms' else '%d' %(int(split_line[1]) + bp_total) if reps_are_chroms_or_chunks == 'chunks' else None
                        f_out_vcf.write('\t'.join(split_line))
    return None



def split_VCF(basefile_vcf, basefile_vcfout, mask_not_just_in = [], mask_not_het_in = [], rename_contig_dict = {}):
    #This splits a simulated VCF into its component chromosomes.
    #I can optionally mask out cases that are
    #i) only not '0|0' for one individual [e.g. H erectus, often 'msp_67'] or
    #ii) heterozygote in one of several individuals [e.g. het in sampled archaic, often 'msp_65','msp_66']

    #rename_contig_dict can be used to rename the contigs.
    
    #First, get the headers
    open_vcf = gzip.open if basefile_vcf[-3:] == '.gz' else open
    open_vcfout = gzip.open if basefile_vcfout[-3:] == '.gz' else open
    headers = []
    with open_vcf(basefile_vcf, 'rb') as f_vcf:
        for line in f_vcf:
            if line[0:1] == '#':
                headers.append(line)
            else:
                break
    #Use the contigs to determine how to split
    contigs = []
    for header in headers:
        if '##contig=' in header:
            contigs.append(header.split('<ID=')[1].split(',')[0])
    #Get IDx of masking individuals
    header_inds = np.array(headers[-1][0:-1].split('\t'))
    num_inds = len(header_inds) - 9
    mask_not_just_in_idx = [np.where(header_inds == i)[0][0] for i in mask_not_just_in]
    
    mask_not_just_in_idx_other = range(9, num_inds + 9)
    for idx in mask_not_just_in_idx:
        mask_not_just_in_idx_other.remove(idx)
    
    mask_not_het_in_idx = [np.where(header_inds == i)[0][0] for i in mask_not_het_in]

    #Read the whole thing in, masking as I go...
    final_line = None
    curr_contig = None
    reading_state = 'reading'
    while reading_state != 'finished':
        lines_read = 0
        vcf_read = []
        with open_vcf(basefile_vcf, 'rb') as f_vcf:
            for line in f_vcf:
                if reading_state == 'writing':
                    break
                if lines_read % 100000 == 1:
                    print("Read in %d SNPs" %(lines_read))
                if lines_read == 1e6: #Work with 1 million lines at a time
                    reading_state = 'finishing'
                    curr_contig = line[:-1].split('\t')[0]
                if reading_state == 'starting':
                    if line == final_line:
                        reading_state = 'reading'
                if reading_state == 'reading' or reading_state == 'finishing':
                    if line[0] == '#':
                        pass
                    else:
                        split_line = line[:-1].split('\t')
                        if reading_state == 'finishing' and split_line[0] != curr_contig:
                            final_line = line
                            reading_state = 'writing'
                            break
                        else:
                            mask = False
                            if len(mask_not_just_in_idx) >= 1:
                                ind_0s, ind_1s = [0,0]
                                for ind_idx in mask_not_just_in_idx_other:
                                    ind_0s += split_line[ind_idx].count('0')
                                    ind_1s += split_line[ind_idx].count('1')
                                not_just_in_0s, not_just_in_1s = [0,0]
                                for ind_idx in mask_not_just_in_idx:
                                    not_just_in_0s += split_line[ind_idx].count('0')
                                    not_just_in_1s += split_line[ind_idx].count('1')
                                if (ind_0s == 0 and not_just_in_0s > 0) or (ind_1s == 0 and not_just_in_1s > 0):
                                    mask = True
                                else:
                                    pass
                                
                            if mask == False and len(mask_not_het_in_idx) >= 1:
                                for ind_idx in mask_not_het_in_idx:
                                    if split_line[ind_idx].count('0') > 0 and split_line[ind_idx].count('1') > 0:
                                        #0|1,1|0,0/1,1/0
                                        mask = True
                                        break
                                    else:
                                        pass
                            
                            if mask == False:
                                vcf_read.append(split_line)
                                lines_read += 1
        
        
        reading_state = 'finished' if reading_state != 'writing' else 'starting'
        vcf_read = np.array(vcf_read)
        contigs_read = np.unique(vcf_read[::,0])
        print("Writing out %d SNPs as %d contigs (%s)" %(lines_read, len(contigs_read), ','.join(contigs_read)))
        for contig in contigs:
            if contig not in rename_contig_dict.keys():
                rename_contig_dict[contig] = contig
            if contig in contigs_read:
                print('writing %s as %s' %(contig, rename_contig_dict[contig]))
                with open_vcfout(basefile_vcfout %(rename_contig_dict[contig]), 'wb') as f_vcfout:
                    for line in headers:
                        if '##contig=' in line:
                            if '<ID=%s,' %(contig) in line:
                                f_vcfout.write(line.replace('<ID=%s,' %(contig), '<ID=%s,' %(rename_contig_dict[contig])))
                        else:
                            f_vcfout.write(line)
                    for line in vcf_read[vcf_read[::,0] == contig]:
                        write_line = copy.copy(line)
                        write_line[0] = rename_contig_dict[contig]
                        f_vcfout.write('\t'.join(write_line) + '\n')
    return None

def combine_mismatch_as_individuals(basefile_mismatch, basefile_out, reps_to_combine, ind_names):
    # This combines all blocks, mismatch or mrca files, splitting them into individuals.
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

def combine_simulation_output(basefile_sims, basefolder_in, basefolder_out, reps_to_combine, ind_names, vcf_reps_are_chroms = False, vcf_reps_are_chunks = False, vcf_reps_are_chunks_insets = None, rename_contig_and_filter_rep_vcfs = False, blocks = False, mismatch = False, mrca = False, sources = []):
    # Simple calling function that combines the data of the VCF files as well as, optionally, blocks, mismatch and mrca files
    # e.g. for combining some H erectus sims try combine_simulation_output(basefile_sims = '/media/sf_Dropbox/Transfer/chunk_analyses/sims/power_sims_herectus/msprime_Malas2016Her.70af60au.heT80000.heN13249.heI02.heTi1353.10Mb.rep_%s', reps_to_combine = [0,1,2,3,4,5,6,7,8,9], ind_names = ['msp_%d' %(i) for i in range(35,65)], blocks = True, mismatch = True, mrca = True, sources = ['deni', 'nean', 'her', 'human'])
    if vcf_reps_are_chroms == True:
        combine_VCF(basefile_vcf = basefolder_in + basefile_sims + '.vcf.gz', basefile_out = basefolder_out + basefile_sims + '.vcf.gz', reps_to_combine = reps_to_combine, reps_are_chroms_or_chunks = 'chroms')
    if vcf_reps_are_chunks == True:
        if type(vcf_reps_are_chunks_insets) != type(None):
            combined = 0
            chromsets = []
            chromset = []
            while combined < len(reps_to_combine):
                chromset.append(reps_to_combine[combined])
                combined += 1
                if len(chromset) == vcf_reps_are_chunks_insets:
                    chromsets.append(chromset)
                    chromset = []
            if len(chromset) > 1:
                chromsets.append(chromset)
            for chrom_idx in range(len(chromsets)):
                combine_VCF(basefile_vcf = basefolder_in + basefile_sims + '.vcf.gz', basefile_out = basefolder_out + basefile_sims %('%schunks_' + str(chrom_idx)) + '.vcf.gz', reps_to_combine = chromsets[chrom_idx], reps_are_chroms_or_chunks = 'chunks')
        else:
            combine_VCF(basefile_vcf = basefolder_in + basefile_sims + '.vcf.gz', basefile_out = basefolder_out + basefile_sims %('%schunks') + '.vcf.gz', reps_to_combine = reps_to_combine, reps_are_chroms_or_chunks = 'chunks')
    if rename_contig_and_filter_rep_vcfs == True:
        for rep in reps_to_combine:
            print("Filtering and renaming rep %s" %(rep))
            raise RuntimeError("This combination method conditions on het status etc in a way that is likely specific for the H erectus case. Revisit!")
            split_VCF(basefile_vcf = basefolder_in + basefile_sims %(rep) + '.vcf.gz', basefile_vcfout = basefolder_out + basefile_sims %('FILTERED_%s') + '.vcf.gz', mask_not_het_in = ['msp_65', 'msp_66'], mask_not_just_in = ['msp_67'], rename_contig_dict = {'1':rep})
    for source in sources:
        if blocks == True:
            combine_mismatch_as_individuals(basefile_mismatch = basefolder_in + basefile_sims + '.blocks.source-%s.BED.gz' %(source), basefile_out = basefolder_out + basefile_sims + '.blocks.source-%s.BED.gz' %(source), reps_to_combine = reps_to_combine, ind_names = ind_names)
        if mrca == True:
            combine_mismatch_as_individuals(basefile_mismatch = basefolder_in + basefile_sims + '.mrca.source-%s.BED.gz' %(source), basefile_out = basefolder_out + basefile_sims + '.mrca.source-%s.BED.gz' %(source), reps_to_combine = reps_to_combine, ind_names = ind_names)
        if mismatch == True:
            combine_mismatch_as_individuals(basefile_mismatch = basefolder_in + basefile_sims + '.mismatch.source-%s.BED.gz' %(source), basefile_out = basefolder_out + basefile_sims + '.mismatch.source-%s.BED.gz' %(source), reps_to_combine = reps_to_combine, ind_names = ind_names)
    return None



