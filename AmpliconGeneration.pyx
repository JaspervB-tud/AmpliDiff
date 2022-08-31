import cython
import numpy as np
from multiprocessing import shared_memory
from Amplicon import *
from math import ceil
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def determine_differences_cy(amplicons, sequences, lineages, ids, int amplicon_threshold, comparison_matrix):
    diffs_per_amp = {a : set() for a in amplicons}

    cdef int seq1, seq2, next_start, i

    for seq1 in range(len(sequences)):
        for seq2 in range(seq1):
            if lineages[seq1] != lineages[seq2]:
                prev_amplicon = (0,0)
                diffs = []
                for next_amplicon in amplicons:
                    diffs = [d for d in diffs if d >= next_amplicon[0]]
                    next_start = max(prev_amplicon[1], next_amplicon[0])
                    for i in range(next_start, next_amplicon[1]):
                        if not comparison_matrix[(sequences[seq1][i], sequences[seq2][i])][0]:
                            diffs.append(i)
                    if len(diffs) > amplicon_threshold:
                        diffs_per_amp[next_amplicon].add((ids[seq2], ids[seq1]))
                    prev_amplicon = next_amplicon
    return diffs_per_amp

@cython.boundscheck(False)
@cython.wraparound(False)
def determine_differences_cy2(amplicons, sequences, lineages, ids, int amplicon_threshold, comparison_matrix, long[:]to_check):
    diffs_per_amp = {a : set() for a in amplicons}

    cdef int seq1, seq2, next_start, i, cur_sum, cur_start, cur_end, last_update
    cdef int[:] cur_diffs_incl
    cdef signed char[:] diffs
    cdef int seq_length = len(sequences[0])

    for seq1 in range(len(sequences)):
        for seq2 in range(seq1):
            if lineages[seq1] != lineages[seq2]:
                cur_sum = 0
                last_update = 0
                cur_diffs_incl = np.zeros((seq_length), dtype=np.int32)
                diffs = np.zeros((seq_length), dtype=np.int8)
                for i in to_check:
                    if not comparison_matrix[(sequences[seq1][i], sequences[seq2][i])][0]:
                        cur_diffs_incl[last_update:i] = cur_sum
                        diffs[i] = 1
                        cur_sum += 1
                        last_update = i
                cur_diffs_incl[last_update:] = cur_sum
                for amplicon in amplicons:
                    cur_start = amplicon[0]
                    cur_end = amplicon[1]-1
                    if cur_diffs_incl[cur_end] - cur_diffs_incl[cur_start] + diffs[cur_start] > amplicon_threshold:
                        diffs_per_amp[amplicon].add((ids[seq2], ids[seq1]))
    return diffs_per_amp

def shm_test(amplicons, sequences, sequence_pairs, str mem_name, int row_index):
    shm = shared_memory.SharedMemory(name=mem_name)
    cdef signed char[:,:,:] np_array = np.ndarray(shape=(len(amplicons), len(sequences), len(sequences)), dtype=np.int8, buffer=shm.buf)

    cdef int seq1, seq2

    for pair in sequence_pairs:
        np_array[row_index][pair[0]][pair[1]] = 1
    print(np.sum(np_array[row_index][:][:]))

def determine_differences_sequencewise_cy(amplicons, sequence_pairs, sequences, lineages, ids, int amplicon_threshold, comparison_matrix, relevant_nucleotides, str mem_name):
    shm = shared_memory.SharedMemory(name=mem_name)
    cdef signed char[:,:,:] diff_tensor = np.ndarray((len(amplicons), len(sequences), len(sequences)), dtype=np.int8, buffer=shm.buf)
    cdef int[:] cur_diffs_cum
    cdef signed char[:] diffs
    cdef int seq1, seq2, cur_sum, last_update, i, cur_start, cur_end, cur_amplicon
    cdef int cur_first = 0
    cdef int seq_length = len(sequences[0])

    for (seq1, seq2) in sequence_pairs:
        if lineages[seq1] != lineages[seq2]:
            cur_sum = 0
            last_update = 0
            cur_diffs_cum = np.zeros((seq_length), dtype=np.int32)
            diffs = np.zeros((seq_length), dtype=np.int8)

            for i in relevant_nucleotides:
                if not comparison_matrix[(sequences[seq1][i], sequences[seq2][i])][0]:
                    cur_diffs_cum[last_update:i] = cur_sum
                    diffs[i] = 1
                    cur_sum += 1
                    last_update = i
            cur_diffs_cum[last_update:] = cur_sum
            cur_amplicon = 0
            for amplicon in amplicons:
                cur_start = amplicon[0]
                cur_end = amplicon[1] - 1
                if cur_diffs_cum[cur_end] - cur_diffs_cum[cur_start] + diffs[cur_start] > amplicon_threshold:
                    diff_tensor[cur_amplicon][seq1][seq2] = 1
                    #diff_tensor[cur_amplicon][seq2][seq1] = 1
                cur_amplicon += 1

def generate_amplicons_sequencewise_cy(amplicons, included_amplicons, str mem_name, int num_sequences, int num_amplicons):
    shm = shared_memory.SharedMemory(name=mem_name)
    diff_tensor = np.ndarray((num_amplicons, num_sequences, num_sequences), dtype=np.int8, buffer=shm.buf)

    res = []
    
    cdef int amplicon_index = 0
    cdef signed char comp = 1
    for amplicon in amplicons:
        cur_diffs = np.where(diff_tensor[included_amplicons[amplicon_index]][:][:] == comp)
        res.append(Amplicon(amplicon[0], amplicon[1]))
        for diff in range(len(cur_diffs[0])):
            res[-1].differences.add((cur_diffs[0][diff],cur_diffs[1][diff]))
        amplicon_index += 1

    return res