import cython
import numpy as np
cimport cython

@cython.wraparound(False)
def generate_amplicons_cy(int[:,:] AMPS, int amplicon_width, int num_amps, signed char[:,:] sequences, int[:,:] sequence_pairs, int total_sequence_pairs, int num_sequences, int[:] ids, signed char[:,:] comparison_matrix, long[:] relevant_nucleotides, int num_relevant, int amplicon_threshold):
    cdef signed char[:,:,:] diff_tensor = np.zeros((num_amps, num_sequences, num_sequences), dtype=np.int8)

    cdef int[:] diffs_cum
    cdef int seq1, seq2, amp, cur_sum, j, cur_index, cur_lb, cur_ub
    cdef long i

    for cur_pair in range(total_sequence_pairs):
        seq1 = sequence_pairs[cur_pair][0]
        seq2 = sequence_pairs[cur_pair][1]
        cur_sum = 0
        diffs_cum = np.zeros((num_relevant), dtype=np.int32)
        cur_index = 0
        for j in range(num_relevant):
            i = relevant_nucleotides[j]
            cur_sum += comparison_matrix[sequences[seq1][i], sequences[seq2][i]]
            diffs_cum[cur_index] = cur_sum
            cur_index += 1
        for amp in range(num_amps):
            cur_lb = relevant_nucleotides[AMPS[amp][1]]
            cur_ub = relevant_nucleotides[AMPS[amp][2]]
            if AMPS[amp][0] <= cur_lb:
                if diffs_cum[AMPS[amp][2]] > amplicon_threshold:
                    #diffs_per_amp[(AMPS[amp][0], AMPS[amp][0] + amplicon_width)].append((ids[seq1],ids[seq2]))
                    diff_tensor[amp][ids[seq1]][ids[seq2]] = 1
            else:
                if diffs_cum[AMPS[amp][2]] - diffs_cum[AMPS[amp][1]] > amplicon_threshold:
                    #diffs_per_amp[(AMPS[amp][0], AMPS[amp][0] + amplicon_width)].append((ids[seq1],ids[seq2]))
                    diff_tensor[amp][ids[seq1]][ids[seq2]] = 1
    return diff_tensor