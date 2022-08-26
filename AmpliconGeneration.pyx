import cython

def move_window_exp_cy(str seq1, str seq2, diffs, prev_amplicon, next_amplicon, comparison_matrix):
    cdef int i

    if prev_amplicon[1] > next_amplicon[0]:
        if len(diffs) > 0:
            while diffs[0] < next_amplicon[0]:
                diffs.pop(0)
                if len(diffs) == 0:
                    break
        for i in range(prev_amplicon[1], next_amplicon[1]):
            if not comparison_matrix[(seq1[i], seq2[i])][0]:
                diffs.append(i)
    else:
        diffs = []
        for i in range(next_amplicon[0], next_amplicon[1]):
            if not comparison_matrix[(seq1[i], seq2[i])][0]:
                diffs.append(i)
    return diffs

def determine_differences_exp_cy(amplicons, sequences, lineages, ids, int amplicon_threshold, comparison_matrix):
    diffs_per_amp = {a : set() for a in amplicons}
    
    cdef int seq1, seq2
    for seq1 in range(len(sequences)):
        for seq2 in range(seq1):
            if lineages[seq1] == lineages[seq2]:
                pass
            else:
                prev_amplicon = (0,0)
                diffs = []
                for next_amplicon in amplicons:
                    diffs = move_window_exp_cy(sequences[seq1], sequences[seq2], diffs, prev_amplicon, next_amplicon, comparison_matrix)
                    if len(diffs) > amplicon_threshold:
                        diffs_per_amp[next_amplicon].add((ids[seq2], ids[seq1]))
                    prev_amplicon = next_amplicon
    return diffs_per_amp