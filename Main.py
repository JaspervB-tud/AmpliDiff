from Generic_Methods import *
from Class_Methods import *

#Can be removed later
import AmpliconGeneration
import time
import numpy as np
import multiprocessing
from math import ceil
from multiprocessing import shared_memory

def sm_test(sequences, amplicons):
    OG_copy = np.zeros((len(amplicons), len(sequences), len(sequences)), dtype=np.int8)

    shm = shared_memory.SharedMemory(create=True, size=(len(amplicons), len(sequences), len(sequences)))
    dst = np.ndarray(OG_copy.shape, dtype=np.int8, buffer=shm.buf)
    dst[:] = OG_copy[:]

    print(dst.nbytes / 10**9)

    sequence_pairs = list(itertools.combinations(sequences,2))
    seqs_partitioned = [sequence_pairs[i:i+ceil(len(sequence_pairs)/4)] for i in range(0, len(sequence_pairs), ceil(len(sequence_pairs)/4))]

    with mp.Pool(4) as pool:
        pool.starmap(AmpliconGeneration.shm_test, zip(itertools.repeat(amplicons), itertools.repeat(sequences), seqs_partitioned, itertools.repeat(shm.name), [1,2,3,4]))

    shm.close()
    shm.unlink()

def translate_to_numeric(sequences, amplicons, relevant_nucleotides, comparison_matrix):
    
    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']
    char_comp = np.zeros((len(chars), len(chars)), dtype=np.int8)
    chars2num = {}
    seqs_num = np.zeros((len(sequences), sequences[0].length), dtype=np.int32)
    amplicons_num = np.zeros((len(amplicons)), dtype=np.int32)
    amplicons_lb = np.zeros((len(amplicons)), dtype=np.int32)
    amplicons_ub = np.zeros((len(amplicons)), dtype=np.int32)
    
    for char_index in range(len(chars)):
        chars2num[chars[char_index]] = char_index
        
    for c1 in range(len(chars)):
        for c2 in range(len(chars)):
            if not comparison_matrix[(chars[c1], chars[c2])][0]:
                char_comp[c1][c2] = 1
                
    for s in range(len(sequences)):
        for i in range(sequences[s].length):
            seqs_num[s][i] = chars2num[sequences[s].sequence[i]]
            
    for a in range(len(amplicons)):
        amplicons_num[a] = amplicons[a][0]
        cur = np.where(relevant_nucleotides < amplicons[a][0])[0]
        if cur.shape[0] > 0:
            amplicons_lb[a] = cur[-1]
        cur = np.where(relevant_nucleotides < amplicons[a][1])[0]
        if cur.shape[0] > 0:
            amplicons_ub[a] = cur[-1]
                
    return chars2num, char_comp, seqs_num, amplicons_num, amplicons_lb, amplicons_ub


if __name__ == '__main__':
    sequences = generate_sequences('/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing', '/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing')
    sequences = sequences[:500]
    sequences, lb, ub, feasible_amplicons, relevant_nucleotides = preprocess_sequences(sequences, 50, amplicon_width=400, misalign_threshold=10)
    
    seqs = [s.sequence for s in sequences]
    lineages = [s.lineage_num for s in sequences]
    ids = [s.id_num for s in sequences]

    amps = list(feasible_amplicons)
    amps.sort(key = lambda x : x[0])
    amps = amps[:10000]

    amplicon_width = 400
    M = generate_opportunistic_matrix()
    c, C, S, A, A_lb, A_ub = translate_to_numeric(sequences, amps, relevant_nucleotides, M)
    
    rel = relevant_nucleotides[relevant_nucleotides < amps[-1][1]]

    runtimes = 0

    st = time.time()
    generate_amplicons_mp_smartest(sequences, 400, M, amplicon_threshold=1, feasible_amplicons=amps, relevant_nucleotides=rel, processors=4)
    print('New runtime: ' + str(time.time() - st))

    st = time.time()
    generate_amplicons_mp_exp_cy(sequences, 400, M, amplicon_threshold=1, feasible_amplicons=amps, processors=4)
    print('Old runtime: ' + str(time.time() - st))

    '''
    same = True
    if len(amps1) != len(amps2):
        same = False
        print('Unequal lengths: ' + str(len(amps1)) + ', ' + str(len(amps2)))
    for i in range(len(amps1)):
        if amps1[i].differences != amps2[i].differences:
            print(str(len(amps1[i].differences)) + ', ' + str(len(amps2[i].differences)))
            same = False
    print(same)

    st = time.time()
    amps1 = AmpliconGeneration.determine_differences_cy(amps, seqs, lineages, ids, 1, M)
    print('Initial method: ' + str(time.time() - st))
    print(len(amps1))
    

    st = time.time()
    amps1 = generate_amplicons_mp_exp_cy(sequences[:500], 400, M, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amps, processors=2)
    print('Split on amplicons: ' + str(time.time() - st))
    
    st = time.time()
    res = []
    amps2 = AmpliconGeneration.generate_amplicons_smarter_cy(A, A_lb, A_ub, amplicon_width, A.shape[0], S, S.shape[0], lineages, C, rel, rel.shape[0], 1)
    for amp in amps2:
        res.append(Amplicon(amp[0], amp[1]))
        res[-1].differences = amps2[amp]
    print('New method: ' + str(time.time() - st))
    print(len(amps2))

    same = True
    print(len(amps1) == len(res))
    for i in range(len(amps1)):
        if amps1[i].differences != res[i].differences:
            same = False
    print(same)


    for a in amps1:
        if amps1[a] != amps2[a]:
            same = False
            print(a)
            print(amps1[a] - amps2[a])
            print('-'*200)
            print(amps2[a] - amps1[a])
            #print(str(len(amps1[a])) + ', ' + str(len(amps2[a])))
            break
    print(same)
    '''
    '''
    n_seqs = 100
    n_amps = 2000
    amplicon_threshold = 1
    comparison_matrix = generate_opportunistic_matrix()

    seqs = [s.sequence for s in sequences]
    amplicons = list(feasible_amplicons)
    amplicons.sort(key = lambda x : x[0])
    lineages = [s.lineage_num for s in sequences]

    ids = [s.id_num for s in sequences]

    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']

    st = time.time()
    #amps1 = AmpliconGeneration.determine_differences_cy(amplicons[:n_amps], seqs[:n_seqs], lineages[:n_seqs], ids[:n_seqs], 1, comparison_matrix)
    amps1 = generate_amplicons_mp_exp_cy(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], processors=4)
    print('Split on amplicons: ' + str(time.time() - st))

    st = time.time()
    amps2 = generate_amplicons_mp_smart(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], relevant_nucleotides=relevant_nucleotides, processors=4)
    print('Split on sequences: ' + str(time.time() - st))
    """
    st = time.time()
    #amps2 = AmpliconGeneration.determine_differences_cy2(amplicons[:n_amps], seqs[:n_seqs], lineages[:n_seqs], ids[:n_seqs], 1, comparison_table, np.where(final_to_check == 1)[0])
    amps2 = generate_amplicons_mp_exp_cy2(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], relevant_nucleotides=relevant_nucleotides, processors=2)
    print(time.time() - st)
    """
    print(len(amps1) == len(amps2))
    same = True
    for i in range(len(amps1)):
        if not amps1[i].differences == amps2[i].differences:
            print(amps1[i].differences - amps2[i].differences)
            print(amps2[i].differences - amps1[i].differences)
            print(len(amps1[i].differences))
            print(len(amps2[i].differences))
            same = False
            break
    print(same)
    '''