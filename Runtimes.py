import Classes
import Scripts
import argparse
import time
import random

def main():
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the greedy amplicon and primer selection algorithm.')
    #Input data
    parser.add_argument('sequences', type=str, help='File location of the (aligned) sequences')
    parser.add_argument('metadata', type=str, help='File location of the metadata for the sequences')
    #Output data
    parser.add_argument('output', type=str, help='Folder where results will be stored')
    #Amplicon parameters
    parser.add_argument('-aw', '--amplicon_width', default=200, type=int, help='Amplicon size')
    parser.add_argument('-at', '--amplicon_threshold', default=0, type=int, help='Number of allowed mismatches in an amplicon')
    parser.add_argument('-mis', '--misalign_threshold', default=5, type=int, help='Number of allowed misalign characters in an amplicon')
    #Primer parameters
    parser.add_argument('-pw', '--primer_width', default=25, type=int, help='Primer size')
    parser.add_argument('-sw', '--search_width', default=50, type=int, help='Search window for finding primers')
    parser.add_argument('-cov', '--coverage', default=1.0, type=float, help='Fraction of sequences that should be covered by both a forward and reverse primer for every amplicon')
    #Greedy parameters
    parser.add_argument('-amps', '--amplicons', default=5, type=int, help='Number of amplicons to find')
    parser.add_argument('-mp', '--multiplex', action='store_true', help='If supplied, will find primers that work in a single batch (multiplex)')
    #General parameters
    parser.add_argument('-c', '--cores', default=1, type=int, help='Number of processor cores to use when using multiprocessing')
    parser.add_argument('-vl', '--variants_location', default=None, type=str, help='File location containing the lineages per variant')
    parser.add_argument('-v', '--variants', default=None, type=str, help='Letters representing the variants to consider')
    parser.add_argument('-n', '--n_sequences', default=10**9, type=int, help='Number of sequences to include')
    parser.add_argument('-s', '--seed', default=1234, type=int, help='Seed to use when randomizing sequences')
    args = parser.parse_args()
    
    #Initialize variables to store information
    runtimes = []
    cur_time = time.time()
    
    #Read sequences
    st = time.time()
    sequences = Scripts.generate_sequences(args.metadata, args.sequences)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent generating sequences: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent generating sequences: ' + str(time.time() - st))
    
    #Generate comparison matrix
    comparison = Scripts.generate_opportunistic_matrix()
    
    #Preprocess sequences
    st = time.time()
    if args.variants_location and args.variants:
        variants = []
        for variant in args.variants:
            if variant == 'a' and 'Alpha' not in variants:
                variants.append('Alpha')
            elif variant == 'b' and 'Beta' not in variants:
                variants.append('Beta')
            elif variant == 'c' and 'Gamma' not in variants:
                variants.append('Gamma')
            elif variant == 'd' and 'Delta' not in variants:
                variants.append('Delta')
            elif variant == 'e' and 'Epsilon' not in variants:
                variants.append('Epsilon')
            elif variant == 'z' and 'Zeta' not in variants:
                variants.append('Zeta')
            elif variant == 'n' and 'Eta' not in variants:
                variants.append('Eta')
            elif variant == 'k' and 'Kappa' not in variants:
                variants.append('Kappa')
            elif variant == 'm' and 'Mu' not in variants:
                variants.append('Mu')
            elif variant == 'o' and 'Omicron' not in variants:
                variants.append('Omicron')
        sequences, lb, ub, feasible_amplicons = Scripts.preprocess_sequences(sequences, args.search_width, variants_location=args.variants_location, variants=variants, amplicon_width=args.amplicon_width, misalign_threshold=args.misalign_threshold)
        #Randomize sequences
        random.seed(args.seed)
        random.shuffle(sequences)
        sequences = sequences[:args.n_sequences]
        with open(args.output + '/sequences_included_' + str(args.seed) + '.txt') as f:
            for sequence in sequences:
                f.write(sequence.id + '\n')
        
        with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
            f.write('Variants considered:\n')
            for variant in variants:
                f.write(variant + '\n')
            f.write('Total sequences = ' + str(len(sequences)) + '\n')
    else:
        sequences, lb, ub, feasible_amplicons = Scripts.preprocess_sequences(sequences, args.search_width, amplicon_width=args.amplicon_width, misalign_threshold=args.misalign_threshold)
        #Randomize sequences
        random.seed(args.seed)
        random.shuffle(sequences)
        sequences = sequences[:args.n_sequences]
        with open(args.output + '/sequences_included_' + str(args.seed) + '.txt') as f:
            for sequence in sequences:
                f.write(sequence.id + '\n')
        
        with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
            f.write('Variants considered:\n')
            for variant in ['Alpha','Beta','Gamma','Delta','Epsilon','Zeta','Eta','Kappa','Mu','Omicron']:
                f.write(variant + '\n')
            f.write('Total sequences = ' + str(len(sequences)) + '\n')
    
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent pre-processing sequences and determining feasible amplicons: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent pre-processing sequences and determining feasible amplicons: ' + str(time.time() - st))
    
    #Generate primer index
    st = time.time()
    PI = Classes.PrimerIndex.generate_index_mp(sequences, args.primer_width, comparison, processors=args.cores)
    PI.remove_redundant()
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent generating primer index and filtering for feasible primers: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent generating primer index and filtering for feasible primers: ' + str(time.time() -st))
    
    #Generate amplicons
    st = time.time()
    amplicons = Scripts.generate_amplicons_mp_exp(sequences, args.amplicon_width, comparison, feasible_amplicons=feasible_amplicons, processors=args.cores, amplicon_threshold=args.amplicon_threshold)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent generating amplicon differentiation ' + str(time.time() - st) + '\n')
        f.write('Total feasible amplicons: ' + str(len(amplicons)) + '\n')
    #runtimes.append('Time spent generating amplicon differentiation: ' + str(time.time() - st))
    
    #Run greedy
    st = time.time()
    logs, amplicons, result_amplicons = Scripts.greedy(sequences, amplicons, args.primer_width, args.search_width, PI, comparison, args.amplicons, args.coverage, 5, logging=True, multiplex=args.multiplex)
    with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
        f.write('Time spent running greedy algorithm: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent running greedy algorithm: ' + str(time.time() - st))
    
    #Run final optimization
    if args.multiplex:
        st = time.time()
        Scripts.check_primer_feasibility(sequences, result_amplicons, PI, optimize=1, coverage=args.coverage)
        with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
            f.write('Time spent doing final primer optimization: ' + str(time.time() - st))
    else:
        st = time.time()
        for amplicon in result_amplicons:
            cur_primers = Scripts.check_primer_feasibility(sequences, [amplicon], PI, optimize=1, coverage=args.coverage)
            with open(args.output + '/runtimes_' + str(args.seed) + '.txt', 'a') as f:
                f.write('Amplicon: ' + str(amplicon.id) + '\n')
                f.write('Forward primers\n')
                for fwd in cur_primers['forward']:
                    f.write(fwd + '\n')
                f.write('Reverse primers\n')
                for rev in cur_primers['reverse']:
                    f.write(rev + '\n')
                
    with open(args.output + '/logfile_' + str(args.seed) + '.txt', 'w') as f:
        for line in logs:
            f.write(line + '\n')
            
    
            
            
            