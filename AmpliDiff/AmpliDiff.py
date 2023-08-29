import argparse
import time
import random
from classless_methods import generate_comparison_matrix
from class_methods import generate_sequences, process_sequences, generate_amplicons, greedy_amplicon_selection
import PrimerIndex

def main():
    parser = argparse.ArgumentParser(description='Run AmpliDiff to find discriminatory amplicons and corresponding primers in input sequences.')
    #High level input data
    parser.add_argument('sequences', type=str, help='Path to aligned sequences fasta file')
    parser.add_argument('metadata', type=str, help='Path to the metadata for the sequences')
    parser.add_argument('-o', '--output', type=str, default='.', help='Path to the folder where output will stored, default is current folder')
    parser.add_argument('--primer_thresholds', type=str, default='./primer_thresholds.csv', help='Path to the primer thresholds file, default is ./primer_thresholds.csv')
    #Amplicon parameters
    parser.add_argument('-aw', '--amplicon_width', type=int, default=200, help='Amplicon size, default is 200')
    parser.add_argument('-mm', '--max_mismatches', type=int, default=0, help='Number of allowed mismatches for amplicon differentiation, default is 0')
    parser.add_argument('-mt', '--max_misalign', type=int, default=20, help='Number of allowed misalign characters in an amplicon, default is 20')
    parser.add_argument('--min_non_align', type=int, default=0, help='Minimum number of nucleotides before the first amplicon and after the final amplicon, default is 0')
    #Primer parameters
    parser.add_argument('-pw', '--primer_width', type=int, default=25, help='Primer size, default is 25')
    parser.add_argument('-sw', '--search_width', type=int, default=50, help='Search window for finding primers, default is 50')
    parser.add_argument('-cov', '--coverage', type=float, default=1.0, help='Minimal required amplifiability, default is 100')
    parser.add_argument('-b', '--beta', type=float, default=0.05, help='Trade-off parameter between primer pairs and differentiability, default is 0.05')
    parser.add_argument('--max_primer_degeneracy', type=int, default=4**5, help='Maximum allowed degeneracy for disambiguating primer candidates, default is 1024')
    parser.add_argument('--gc_lb', type=float, default=0.4, help='Minimum required GC content in primers, default is 0.4')
    parser.add_argument('--gc_ub', type=float, default=0.6, help='Maximum allowed GC content in primers, default is 0.6')
    parser.add_argument('--melting_lb', type=float, default=55., help='Minimum required primer melting temperature in degrees Celsius, default is 55')
    parser.add_argument('--melting_ub', type=float, default=75., help='Maximum allowed primer melting temperature in degrees Celsius, default is 75')
    parser.add_argument('--end_at_threshold', type=int, default=2, help="Maximium allowed A/T nucleotides in final 3 nucleotides, default is 2")
    parser.add_argument('--end_gc_threshold', type=int, default=3, help="Maximum allowed C/G nucleotides in final 5 nucleotides, default is 3")
    parser.add_argument('--monorun_threshold', type=int, default=3, help='Maximum allowed number of single nucleotide runs, default is 3')
    parser.add_argument('--duorun_threshold', type=int, default=3, help="Maximum allowed number of double nucleotide runs, default is 3")
    parser.add_argument('--mfe_threshold', type=float, default=-5., help='Minimum required MFE for determining hairpin formation risk, default is -5')
    parser.add_argument('--self_complementarity_threshold', type=int, default=10, help='Maximum self and primer-primer complementarity in "worst" alignment, default is 10')
    parser.add_argument('--max_temperature_difference', type=float, default=5, help='Maximum difference between minimum and maximum primer melting temperatures, default is 5')
    #Greedy algorithm parameters
    parser.add_argument('-amps', '--num_amplicons', type=int, default=10, help='Number of amplicons to find, default is 10')
    #Sequence parameters
    parser.add_argument('-n', '--num_sequences', type=int, default=-1, help='Number of sequences to include, default is -1 (no maximum)')
    parser.add_argument('--min_characters', type=int, default=-1, help='Minimum number of characters in a sequence, default is -1 (no minimum)')
    parser.add_argument('--max_degeneracy', type=int, default=-1, help='Maximum number of degenerate degeneracy of a sequence, default is -1 (no maximum)')
    parser.add_argument('--max_n', type=int, default=-1, help='Maximum number of n characters in a sequence, default is -1 (no maximum)')
    #System parameters
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of processing cores to use, default is 1 (no multiprocessing)')
    parser.add_argument('-sd', '--seed', type=int, default=0, help='Seed that is used to determine selection of sequences if more than allowed, default is 0')
    
    #Parse arguments
    args = parser.parse_args()
    comparison_matrix = generate_comparison_matrix()

    #Read sequences
    st = time.time()
    print('Reading sequences')
    sequences = generate_sequences(args.sequences, args.metadata, min_characters=args.min_characters,
                                    max_degeneracy=args.max_degeneracy, max_n=args.max_n)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'w+') as f:
        f.write('Time spent generating sequences: ' + str(time.time() - st) + '\n')
    print('Done reading sequences')

    #Randomly pick sequences up until maximum number of sequences to include
    st = time.time()
    print('Randomly selecting up to ' + str(args.num_sequences) + ' sequences with seed=' + str(args.seed))
    random.seed(args.seed)
    random.shuffle(sequences)
    sequences = sequences[:args.num_sequences]
    #Assign new numerical ids to sequences
    for i in range(len(sequences)):
        sequences[i].id_num = i
    print('Done selecting sequences')
    #Save sequence ids of sequences that were included for reproducibility
    with open(args.output + '/sequences_included_' + str(args.seed) + '.txt', 'w+') as f:
        for sequence in sequences:
            f.write(sequence.id + '\n')
    
    #Process sequences
    print('Processing sequences')
    sequences, lb, ub, feasible_amplicons, relevant_nucleotides = process_sequences(sequences, min_non_align=args.min_non_align, 
                                                                                    amplicon_width=args.amplicon_width, max_misalign=args.max_misalign)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent processing sequences and determining feasible amplicons: ' + str(time.time() - st) + ', number of feasible amplicons: ' + str(len(feasible_amplicons)) + '\n')
    print('Done processing sequences')

    #Generate PrimerIndex
    st = time.time()
    #Read thresholds for primers, prioritizing commandline arguments over primer_thresholds.csv
    default_thresholds = {
            'gc_lb'                             : 0.4,
            'gc_ub'                             : 0.6,
            'melting_lb'                        : 55.,
            'melting_ub'                        : 75.,
            'end_at_threshold'                  : 2,
            'end_gc_threshold'                  : 3,
            'monorun_threshold'                 : 3,
            'duorun_threshold'                  : 3,
            'mfe_threshold'                     : -5.,
            'self_complementarity_threshold'    : 10
            }
    thresholds = {}
    try:
        with open(args.primer_thresholds, 'r') as f:
            for threshold in f:
                threshold = threshold.strip().split(';')
                if threshold[0] in ['gc_lb', 'gc_ub', 'melting_lb', 'melting_ub', 'mfe_threshold']:
                    thresholds[threshold[0]] = float(threshold[1])
                elif threshold[0] in ['end_at_threshold', 'end_gc_threshold', 'monorun_threshold', 'duorun_threshold', 'self_complementarity_threshold']:
                    thresholds[threshold[0]] = int(threshold[0])
    except Exception as e:
        print("Couldn't read primer_thresholds.csv, using default thresholds or those supplied as command line arguments")
        thresholds = default_thresholds
    if thresholds['gc_lb'] == default_thresholds['gc_lb']:
        thresholds['gc_lb'] = args.gc_lb
    if thresholds['gc_ub'] == default_thresholds['gc_ub']:
        thresholds['gc_ub'] = args.gc_ub
    if thresholds['melting_lb'] == default_thresholds['melting_lb']:
        thresholds['melting_lb'] = args.melting_lb
    if thresholds['melting_ub'] == default_thresholds['melting_ub']:
        thresholds['melting_ub'] = args.melting_ub
    if thresholds['end_at_threshold'] == default_thresholds['end_at_threshold']:
        thresholds['end_at_threshold'] = args.end_at_threshold
    if thresholds['end_gc_threshold'] == default_thresholds['end_gc_threshold']:
        thresholds['end_gc_threshold'] = args.end_gc_threshold
    if thresholds['monorun_threshold'] == default_thresholds['monorun_threshold']:
        thresholds['monorun_threshold'] = args.monorun_threshold
    if thresholds['duorun_threshold'] == default_thresholds['duorun_threshold']:
        thresholds['duorun_threshold'] = args.duorun_threshold
    if thresholds['mfe_threshold'] == default_thresholds['mfe_threshold']:
        thresholds['mfe_threshold'] = args.mfe_threshold
    if thresholds['self_complementarity_threshold'] == default_thresholds['self_complementarity_threshold']:
        thresholds['self_complementarity_threshold'] = args.self_complementarity_threshold
    PrimerIndex.PrimerIndex.set_thresholds(thresholds)
    print('Generating primer index')
    primer_index = PrimerIndex.PrimerIndex.generate_index_mp(sequences, args.primer_width, comparison_matrix, max_degeneracy=args.max_primer_degeneracy, processors=args.cores)
    primer_index.remove_redundant()
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent generating primer index and filtering for feasible primers: ' + str(time.time() - st) + '\n')
    print('Done generating primer index')

    #Generate amplicons
    st = time.time()
    print('Determining amplicon differentiabilities')
    amplicons, differences_per_amplicon = generate_amplicons(sequences, args.amplicon_width, comparison_matrix, 
                        max_mismatch=args.max_mismatches, feasible_amplicons=feasible_amplicons, relevant_nucleotides=relevant_nucleotides)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent determining amplicon differentiabilities: ' + str(time.time() - st) + '\n')
    print('Done determining amplicon differentiabilities')

    #Run greedy algorithm
    st = time.time()
    print('Running greedy algorithm')
    logs, result_amplicons, result_primers = greedy_amplicon_selection(sequences, amplicons, differences_per_amplicon, 
                                args.primer_width, args.search_width, primer_index, comparison_matrix,
                                args.num_amplicons, args.coverage, args.max_temperature_difference, logging=True,
                                beta=args.beta, output_file=args.output+'/primers_'+str(args.seed)+'.fasta')
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent running greedy algorithm: ' + str(time.time() - st) + '\n')
    print('Done running greedy algorithm')

    #Finalizing output
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        n_fwd = 0
        n_rev = 0
        cur = 0
        for primer_set in result_primers:
            f.write('Amplicon: ' + str(result_amplicons[cur].id) + '\n')
            f.write('Forward primers\n')
            for primer in primer_set['forward']:
                f.write('>F' + str(n_fwd) + '\n')
                f.write(primer + '\n')
                n_fwd += 1
            for primer in primer_set['reverse']:
                f.write('>R' + str(n_rev) + '\n')
                f.write(primer + '\n')
                n_rev += 1
            cur += 1
    #Write logfile
    with open(args.output + '/logfile_' + str(args.seed) + '.txt', 'w') as f:
        for line in logs:
            f.write(line + '\n')

if __name__ == '__main__':
    main()