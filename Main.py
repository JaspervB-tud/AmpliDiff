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


if __name__ == '__main__':
    sequences = generate_sequences('/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing', '/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing')
    sequences = sequences[:1000]
    sequences, lb, ub, feasible_amplicons, relevant_nucleotides = preprocess_sequences(sequences, 50, amplicon_width=400, misalign_threshold=10)

    n_seqs = 1000
    n_amps = 10000
    amplicon_threshold = 1
    comparison_matrix = generate_opportunistic_matrix()

    seqs = [s.sequence for s in sequences]
    amplicons = list(feasible_amplicons)
    amplicons.sort(key = lambda x : x[0])
    lineages = [s.lineage_num for s in sequences]

    ids = [s.id_num for s in sequences]

    

    st = time.time()
    #amps1 = AmpliconGeneration.determine_differences_cy(amplicons[:n_amps], seqs[:n_seqs], lineages[:n_seqs], ids[:n_seqs], 1, comparison_matrix)
    amps1 = generate_amplicons_mp_exp_cy(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], processors=5)
    print('Split on amplicons: ' + str(time.time() - st))

    del amps1

    st = time.time()
    amps2 = generate_amplicons_mp_sequences(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], relevant_nucleotides=relevant_nucleotides, processors=5)
    print('Split on sequences: ' + str(time.time() - st))
    """
    st = time.time()
    #amps2 = AmpliconGeneration.determine_differences_cy2(amplicons[:n_amps], seqs[:n_seqs], lineages[:n_seqs], ids[:n_seqs], 1, comparison_table, np.where(final_to_check == 1)[0])
    amps2 = generate_amplicons_mp_exp_cy2(sequences[:n_seqs], 200, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=amplicons[:n_amps], relevant_nucleotides=relevant_nucleotides, processors=2)
    print(time.time() - st)

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
    """