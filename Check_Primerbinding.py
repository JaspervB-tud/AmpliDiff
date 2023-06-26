import multiprocessing as mp
import itertools
from Fetch_Amplicons import generate_sequences, read_logfile, read_primerfile
from Generic_Methods import reverse_complement, generate_opportunistic_matrix
import time
import os
import argparse

def calculate_primer_binding(sequence, ident, amplicons, opportunistic_matrix):
    return (ident, locate_amplicons(sequence, amplicons, opportunistic_matrix, max_degen=2))

def pull_genomes_from_folder(folder):
    sequences = []
    for file in os.listdir(folder + '/Sequences'):
        cur_id = file.split('.')[0]
        if len(cur_id) > 0:
            cur_sequences = generate_sequences(folder + '/Sequences/' + cur_id + '.sequences.fasta', folder + '/Metadata/' + cur_id + '.metadata.tsv')
            sequences.extend(cur_sequences)
    return sequences

def calculate_bindings_per_lineage(sequences, bindings):
    amplicon_bindings = {}
    primer_bindings = {}
    for binding in bindings:
        cur_sequence = sequences[sequences.index(binding[0])]
        if cur_sequence.lineage not in amplicon_bindings:
            amplicon_bindings[cur_sequence.lineage] = {}
            primer_bindings[cur_sequence.lineage] = {}
            for amplicon in binding[1][0]:
                amplicon_bindings[cur_sequence.lineage][amplicon] = 0
                primer_bindings[cur_sequence.lineage][amplicon] = {'forward' : {}, 'reverse' : {}}
                for forward in binding[1][1][amplicon]['forward']:
                    primer_bindings[cur_sequence.lineage][amplicon]['forward'][forward] = 0
                for reverse in binding[1][1][amplicon]['reverse']:
                    primer_bindings[cur_sequence.lineage][amplicon]['reverse'][reverse] = 0
        for amplicon in binding[1][0]:
            if binding[1][0][amplicon]:
                amplicon_bindings[cur_sequence.lineage][amplicon] += 1
            for forward in binding[1][1][amplicon]['forward']:
                if len(binding[1][1][amplicon]['forward'][forward]) > 0:
                    primer_bindings[cur_sequence.lineage][amplicon]['forward'][forward] += 1
            for reverse in binding[1][1][amplicon]['reverse']:
                if len(binding[1][1][amplicon]['reverse'][reverse]) > 0:
                    primer_bindings[cur_sequence.lineage][amplicon]['reverse'][reverse] += 1
    return amplicon_bindings, primer_bindings
            
        
    
def locate_primers(sequence, primerlist, comparison_matrix, max_degen=10):
    '''
    Function that determines the primer binding sites of given primers in the given sequence by considering every exact binding location while
    allowing for degenerate nucleotides (i.e. if the sequence contains an N, any nucleotide will bind here).

    Parameters
    ----------
    sequence : Sequence
        Sequence object for which to determine primer binding sites.
    primerlist : dict{ 'forward':[], 'reverse':[] }
        Dictionary containing a forward and reverse key along with primers (strings) for which binding sites have to be determined.
    comparison_matrix : dict{ char : [] }
        Dictionary which dictates which nucleotides are considered equal.
    max_degen : int, optional
        Number of degenerate nucleotides allowed when checking if a primer would bind. The default is 10 (log2(4**5)).

    Returns
    -------
    hits_fwd : dict{ str : set() }
        Dictionary with the binding sites for every forward primer in $primerlist.
    hits_rev : dict{ str : set() }
        Dictionary with the binding sites for every reverse primer in $primerlist.
    '''
    sequence_reverse = reverse_complement(sequence, rev=True)
    num_degen = len(sequence)-sequence.count('a')-sequence.count('c')-sequence.count('g')-sequence.count('t')
    
    def find_hits(sequence, primers, num_degen, comparison_matrix, reverse=False, max_degen=10):
        '''
        Function that performs degenerate elastic string matching by doing k-mer based matching.

        Parameters
        ----------
        sequence : Sequence
            Sequence object for which to determine primer binding sites.
        primers : list[ str ]
            List with primers for which binding sites have to be determined.
        num_degen : int
            Number of degenerate nucleotides in the sequence. This guides the matching approach.
        comparison_matrix : dict{ char : [] }
            Dictionary which dictates which nucleotides are considered equal.
        reverse : bool, optional
            Boolean which should be set to true if the primers are based on the reverse complement of the sequence. The default is False.
        max_degen : int
            Number of degenerate nucleotides allowed in a primer binding site. The default is 10.

        Returns
        -------
        hits : dict{ str : set() }
            Dictionary with the binding sites for every primer in $primers.

        '''
        hits = {primer: set() for primer in primers}
        equivalent_characters = {'a' : ['a','r','m','w','d','h','v','n'],
                                 'c' : ['c','y','m','s','b','h','v','n'],
                                 'g' : ['g','r','k','s','b','d','v','n'],
                                 't' : ['t','y','k','w','b','d','h','n']}
        sequence_length = len(sequence)
        
        #Behaviour when there are NO degenerate nucleotides
        if num_degen == 0:
            for primer in primers:
                cont = True
                cur_occurrence = -1
                primer_length = len(primer)
                while cont:
                    try:
                        cur_occurrence = sequence.index(primer, cur_occurrence+1)
                        if not reverse:
                            hits[primer].add(cur_occurrence)
                        else:
                            hits[primer].add(sequence_length - cur_occurrence - primer_length)
                    except:
                        cont = False
        #If there are degenerate nucleotides perform a simple version of Boyer-Moore
        else:
            for primer in primers:
                cur_char_index = 0
                primer_length = len(primer)
                while cur_char_index <= sequence_length - primer_length: #iterate over sequence to find occurrence of current primer
                    stepsize = (0, False) #track where we should check for occurrence next
                    cur_degen_chars = 0 #track number of degenerate characters
                    match = True
                    for i in range(primer_length):
                        if sequence[cur_char_index + i] not in ['a','c','g','t']:
                            cur_degen_chars += 1
                        #If current character in sequence is equal to first primer character then store it as next starting point
                        if sequence[cur_char_index+i] in equivalent_characters[primer[0]] and not stepsize[1]:
                            stepsize = (max(1, i), True)
                        #If character does not match then break for-loop and start at next index
                        if not comparison_matrix[(primer[i], sequence[cur_char_index + i])][0]:
                            match = False
                            break
                    if match and cur_degen_chars <= max_degen:
                        if not reverse:
                            hits[primer].add(cur_char_index)
                        else:
                            hits[primer].add(sequence_length - cur_char_index + primer_length)
                    if stepsize[1]:
                        cur_char_index += stepsize[0]
                    else:
                        cur_char_index += i + 1
                        
        return hits
                        
    hits_fwd = find_hits(sequence, primerlist['forward'], num_degen, comparison_matrix, reverse=False)
    hits_rev = find_hits(sequence_reverse, primerlist['reverse'], num_degen, comparison_matrix, reverse=True)
    
    return hits_fwd, hits_rev    
             
def locate_amplicons(sequence, amplicons, comparison_matrix, primer_length=25, max_degen=10):
    '''
    Function that locates the realized amplicons based on the amplicons and corresponding primers in $amplicons in $sequence.
    Note that if the sequence has a stretch of $primer_length degenerate nucleotides, the results can potentially be uninterpretable if $max_degen
    is not configured appropriately.

    Parameters
    ----------
    sequence : Sequence
        Sequence object for which to determine primer binding sites.
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the form of a tuple (start, end), and a dictionary containing the keys forward and reverse which
        contain for every amplicon the corresponding forward and reverse primers in a list as values.
    comparison_matrix : dict{ char : [] }
        Dictionary which dictates which nucleotides are considered equal.
    primer_length : int, optional
        Length of the primer sequence. Theoretically this could be omitted, but since AmpliVar generates primers of fixed length
        it made sense to just include it as a parameter. The default is 25.
    max_degen : int
        Number of degenerate nucleotides in the sequence. This guides the matching approach. The default is 10 (log2(4**5)).

    Returns
    -------
    binding_sites : dict{ (int,int) : (bool) }
        Dictionary with AmpliVar based amplicons as keys, and True if amplicon is amplifiable, False otherwise
    primer_bindings : dict{ (int,int) : dict{ } }
        Dictionary with AmpliDiff based amplicons as keys and a dictionary containing for every corresponding primer the location where it binds on the genome
    '''
    binding_sites = {amplicon[0]: None for amplicon in amplicons}
    primer_bindings = {amplicon[0] : None for amplicon in amplicons}
    #Iterate over amplicons
    for amplicon in amplicons:
        fwd_hits, rev_hits = locate_primers(sequence, amplicon[1], comparison_matrix)
        primer_bindings[amplicon[0]] = {'forward': fwd_hits, 'reverse': rev_hits}
        amplified = (0, 10**16, False)
        
        fwd_indices = set()
        rev_indices = set()
        for fwd in fwd_hits:
            fwd_indices = fwd_indices.union(fwd_hits[fwd])
        for rev in rev_hits:
            rev_indices = rev_indices.union(rev_hits[rev])
        fwd_indices = list(fwd_indices)
        rev_indices = list(rev_indices)
        #print(amplicon[0], fwd_indices, rev_indices)
        
        for fwd, rev in itertools.product(fwd_indices, rev_indices):
            #print(fwd, rev)
            if rev - fwd >= 0 and rev - fwd  < amplified[1] - amplified[0] - primer_length:
                #OLD: Amplicon including primers
                #amplified = (fwd, rev+primer_length, True)
                #NEW: Amplicon excluding primers -> reads and refs will not contain the primers
                amplified = (fwd+primer_length-1, rev, True)
        if amplified[2]:
            binding_sites[amplicon[0]] = True
            
    return (binding_sites, primer_bindings)

def main():
    parser = argparse.ArgumentParser(description='Given an amplicon file, primer file and a folder with genomes determine how often every primer and every amplicon can be amplified in the supplied genomes.')
    #Input data
    parser.add_argument('sequence_folder', type=str, help='Folder containing the genomes (and metadata) for which to check primer bindings and amplicon amplifiability')
    parser.add_argument('logfile', type=str, help='File with the logs for the amplicons')
    parser.add_argument('primerfile', type=str, help='File with the primers corresponding to the amplicons in the logfile')
    #Output data
    parser.add_argument('output', type=str, help='Folder where results will be stored')
    #Additional parameters
    parser.add_argument('-c', '--cores', type=int, help='Number of cores to use')
    
    args = parser.parse_args()
    
    #"Pre-processing"
    M = generate_opportunistic_matrix()
    sequences = pull_genomes_from_folder(args.sequence_folder)
    sequences_raw = []
    ids = []
    genomes_per_lineage = {}
    for sequence in sequences:
        sequences_raw.append(sequence.sequence_raw)
        ids.append(sequence.id)
        if sequence.lineage not in genomes_per_lineage:
            genomes_per_lineage[sequence.lineage] = 0
        genomes_per_lineage[sequence.lineage] += 1
    #Pull amplicons and primers
    amplicons = read_logfile(args.logfile)
    amplicons, primers = read_primerfile(args.primerfile, amplicons)
    #Determine primerbindings
    if args.cores > 1:
        with mp.Pool(args.cores) as pool:
            all_bindings = pool.starmap(calculate_primer_binding, zip(sequences_raw, ids, itertools.repeat(amplicons), itertools.repeat(M)))
    else:
        all_bindings = [calculate_primer_binding(sequence.sequence_raw, sequence.id, amplicons, M) for sequence in sequences]
    #Calculate global stats
    global_amplicon_bindings = {amplicon[0] : 0 for amplicon in amplicons} # stores for every amplicon in how many genomes it is amplifiable
    global_primer_bindings = {amplicon[0] : {'forward' : {}, 'reverse' : {}} for amplicon in amplicons}
    for binding in all_bindings:
        for amplicon in binding[1][0]:
            if binding[1][0][amplicon]:
                global_amplicon_bindings[amplicon] += 1
            for forward in binding[1][1][amplicon]['forward']:
                if forward not in global_primer_bindings[amplicon]['forward']:
                    global_primer_bindings[amplicon]['forward'][forward] = 0
                if len(binding[1][1][amplicon]['forward'][forward]) == 1:
                    global_primer_bindings[amplicon]['forward'][forward] += 1
            for reverse in binding[1][1][amplicon]['reverse']:
                if reverse not in global_primer_bindings[amplicon]['reverse']:
                    global_primer_bindings[amplicon]['reverse'][reverse] = 0
                if len(binding[1][1][amplicon]['reverse'][reverse]) == 1:
                    global_primer_bindings[amplicon]['reverse'][reverse] += 1
    #Calculate stats per lineage
    lineage_amplicon_bindings, lineage_primer_bindings = calculate_bindings_per_lineage(sequences, all_bindings)
    #Start printing to output files
    with open(args.output + '/global_amplicon_stats.csv', 'w') as f:
        f.write('Amplicon;amplifiability;percentage\n')
        for amplicon in amplicons:
            f.write(str(amplicon[0]) + ';' + str(global_amplicon_bindings[amplicon[0]]) + ';' + str(global_amplicon_bindings[amplicon[0]]/len(sequences)) + '\n')
            print(global_amplicon_bindings[amplicon[0]]/len(sequences))
        f.write('Total sequences;' + str(len(sequences)))
    with open(args.output + '/global_primer_stats.csv', 'w') as f:
        f.write('Amplicon/primer;genome binding;percentage\n')
        for amplicon in amplicons:
            f.write('Amplicon:' + str(amplicon[0]) + '\n')
            f.write('forward primers\n')
            for forward in amplicon[1]['forward']:
                f.write(forward + ';' + str(global_primer_bindings[amplicon[0]]['forward'][forward]) + ';' + str(global_primer_bindings[amplicon[0]]['forward'][forward]/len(sequences)) + '\n')
            f.write('reverse primers\n')
            for reverse in amplicon[1]['reverse']:
                f.write(reverse + ';' + str(global_primer_bindings[amplicon[0]]['reverse'][reverse]) + ';' + str(global_primer_bindings[amplicon[0]]['reverse'][reverse]/len(sequences)) + '\n')
        f.write('Total sequences;' + str(len(sequences)))
        
if __name__ == '__main__':
    t = time.time()
    main()
    print(time.time() - t)