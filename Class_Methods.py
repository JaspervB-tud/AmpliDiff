from Bio import AlignIO
import csv
from math import ceil
import multiprocessing as mp
from multiprocessing import shared_memory
import itertools
import time
import gurobipy as gp
from gurobipy import GRB

from Sequence import *
from Amplicon import *
from Generic_Methods import equivalent_characters

import AmpliconGeneration

def generate_sequences(in_folder, out_folder):
    '''
    Function that reads aligned sequences from a "sequences_aligned.fasta" fasta file and saves them as Sequence objects with metadata from "metadata.tsv".

    Parameters
    ----------
    in_folder : str
        Location of the metadata.
    out_folder : str
        Location of the aligned sequences.

    Returns
    -------
    sequences : list[Sequence]
        List of sequences contained in out_folder/sequences_aligned.fasta.

    '''
    sequences_temp = {}
    to_delete = []
    #Read sequences from aligned sequences location which should be in out_folder <- this may be a bit counter-intuitive and could be improved
    aligned_sequence_objects = AlignIO.read(open(out_folder + '/sequences_aligned.fasta'), 'fasta')
    for sequence in aligned_sequence_objects:
        sequences_temp[sequence.id.split('|')[0]] = str(sequence.seq.lower())
        if len(sequence.seq.replace('-','')) < 29000: #Require sequences of at least length 29k
            to_delete.append(sequence.id.split('|')[0]) #store sequences that are probably incorrectly included in the case of SARS-CoV-2        
    #Read metadata, unless impossible in which case we assign every sequence its own "lineage"
    skip = -1
    try:
        sequences = []
        for meta in csv.reader(open(in_folder + '/metadata.tsv'), delimiter='\t'):
            if skip == -1: #First line is the header and shows which column contains the lineage information
                for cur_meta in range(len(meta)):
                    if 'lineage' in meta[cur_meta].lower():
                        skip = cur_meta
                        break
            else:
                #meta[0] always contains the id
                if meta[0] not in to_delete:
                    sequences.append(Sequence(sequences_temp[meta[0].replace(' ', '')], meta[0], lineage=meta[skip]))
                
    except:
        print('Unable to read metadata from file')
        sequences = []
        i = 0
        for identifier in sequences_temp:
            sequences.append(Sequence(sequences_temp[identifier], identifier, lineage=str(i)))
            i += 1
    
    return sequences

def check_variant(sequence, lineages_per_variant, variants):
    '''
    Function that determines whether the sequence is part of one of the variants in $variants.

    Parameters
    ----------
    sequence : Sequence
        Sequence for which we want to know whether it is part of any of the variants of interest (not VOI as defined by WHO).
    lineages_per_variant : dict[ str ]
        Dictionary that contains the lineages corresponding to every variant.
    variants : list[ str ]
        List with variants of interest (not VOI as defined by WHO).

    Returns
    -------
    str
        Variant that the sequence is part of, or None if it is not part of any of the variants in $variants.

    '''
    for variant in variants:
        try:
            for lineage in lineages_per_variant[variant]:
                if lineage in sequence.lineage:
                    return variant
        except:
            continue
    return None

def preprocess_sequences(sequences, min_non_align, variants_location=None, variants=[], amplicon_width=0, misalign_threshold=-1, lineages_location=None, min_sequences_threshold=0):
    '''
    Function that preprocesses the sequences by excluding sequences not part of any of the variants in $variants and
    calculating the lower- and upperbound such that that every sequence contains at least &min_non_align nucleotides before and after
    $lb and $ub respectively. Additionally also determines feasible amplicons given amplicon width and misalignment character threshold, and 
    nucleotides to consider when differentiating sequences.

    Parameters
    ----------
    sequences : list[ Sequence ]
        List with sequences to filter.
    min_non_align : int
        Number of nucleotides to include before (after) $lb ($ub).
    variants_location : str, optional
        Location of the variants.txt file that explains which lineages are part of which variant. The default is None in which case every variant is considered.
    variants : list[ str ], optional
        List of the variants of interest (not VOI as defined by WHO). The default is [] in which case every variant is considered.
    amplicon_width : int, optional
        Size of the amplicons, if we want to determine their feasibility a priori. The default is 0 in which case feasibility of amplicons is not checked.
    misalign_threshold : int, optional
        Number of allowed misalign characters in an amplicon. The default is -1 in which case feasibility of amplicons is not checked.
    lineages_location : str, optional
        Location of the lineages.txt file which contains the number of sequences per lineage. The default is None in which case every lineage is considered.
    min_sequences_threshold : int, optional
        Absolute number of sequences that have to belong to a lineage in order to consider it. The default is 0 in which case every lineage is considered.

    Returns
    -------
    sequences : list[Sequence]
        List of sequences that belong the specified variants.
    lb : int
        Lowerbound such that every sequence has at least $min_non_align nucleotides before it.
    ub : int
        Upperbound such that every sequence has at least $min_non_align nucleotides after it.
    feasible_amplicons : set{ (int,int) }
        Set of feasible amplicons in case amplicon feasibility should be checked.
    relevant_nucleotides : np.array
        Numpy array with the indices of nucleotides that are potentially different between pairs of sequences

    '''
    def determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold):
        '''
        Function that determines feasible amplicons based on the number of misalignment characters.

        Parameters
        ----------
        sequence : Sequence
            Sequence object for which the amplicons should be checked.
        lb : int
            Start index (inclusive) of the first amplicon.
        ub : int
            End index (exclusive) of the final amplicon.
        amplicon_width : int
            Width of the amplicons to check.
        misalign_threshold : int
            Maximum allowed number of misalignment characters.

        Returns
        -------
        feasible_amplicons : set{ (int,int) }
            Set containing the amplicons (start,end) which do not contain too many misalignment characters.

        '''
        feasible_amplicons = set()
        misalign_indices = []

        #Determine misalign indices in initial amplicon minus the final index
        for i in range(lb, lb + amplicon_width - 1):
            if sequence.sequence[i] == '-':
                misalign_indices.append(i)
        #Iterate over all amplicons
        for i in range(lb + amplicon_width - 1, ub):
            #Check if next character is misalign character
            if sequence.sequence[i] == '-':
                misalign_indices.append(i)
            #Check if the current amplicon has too many mis-align characters
            if len(misalign_indices) <= misalign_threshold:
                feasible_amplicons.add((i - amplicon_width + 1, i + 1))
            #Check if the first index in list should be removed
            try:
                if misalign_indices[0] == i - amplicon_width + 1:
                    misalign_indices.pop(0)
            except:
                continue
        return feasible_amplicons
                
    def filter_on_occurrences(sequences, lineages_location, min_sequences_threshold):
        '''
        Function that filters the given sequences based on how often their lineage occurs in the full database.

        Parameters
        ----------
        sequences : list[ Sequence ] 
            List of sequences.
        lineages_location : str
            File location of the lineages.txt file containing occurrences per lineage in the full database.
        min_sequences_threshold : float
            Minimum occurrence of a lineage to be considered.

        Returns
        -------
        list[ Sequence ]
            List of sequences that belong to a lineage that meets the minimum occurrence threshold.

        '''
        filtered_sequences = []
        sequences_per_lineage = {}
        total_sequences = 0
        
        for line in csv.reader(open(lineages_location + '/lineages.txt'), delimiter='\t'):
            try:
                cur_line = line[0].split()
                #Count the number of sequences of each lineage in database
                sequences_per_lineage[cur_line[0]] = int(cur_line[1])
                total_sequences += int(cur_line[1])
            except:
                print('Line: "' + str(line) + '" cannot be parsed and hence is not considered!')
                continue
        print('Total sequences: %d\n' % total_sequences)
        #Check which sequences belong to a lineage that exceeds the minimum threshold
        index = 0 #<- sequence index
        for sequence in sequences:
            #Try-Except in case something weird happens to lineages.txt
            try:
                if sequences_per_lineage[sequence.lineage] >= min_sequences_threshold:
                    filtered_sequences.append(index)
            except:
                print(sequence.id + ' : ' + sequence.lineage + ' does not occur in lineages.txt and is not considered!')
            index += 1
        return [sequences[i] for i in filtered_sequences]

    def filter_on_variants(sequences, variants_location, variants):
        '''
        Function that filters the given sequences based on whether they are part of a lineage that belongs to a variant in $variants.

        Parameters
        ----------
        sequences : list[ Sequence ]
            List of sequences.
        variants_location : str
            File location of the variants.txt file containing the (super-)lineages per variant.
        variants : list[ str ]
            List of variants to include.

        Returns
        -------
        list[ Sequence ]
            List of sequences that belong to a lineage that is part of the variants of interest (not VOI).

        '''
        filtered_sequences = []
        lineages_per_variant = {}
        
        for variant in csv.reader(open(variants_location + '/variants.txt'), delimiter='\t'):
            try:
                lineages_per_variant[variant[0]] = []
                for lineage in variant[1:]:
                    if lineage != '':
                        lineages_per_variant[variant[0]].append(lineage)
            except:
                print('Line "' + str(variant) + '" cannot be parsed and hence is not considered!')
                continue
        #Check which sequences belong to a lineage part of the variants to be considered
        index = 0 #<- sequence index
        for sequence in sequences:
            cur_variant = check_variant(sequence, lineages_per_variant, variants) #check if lineage is part of relevant variants
            if cur_variant:
                print( 'Sequence %s with lineage %s is part of the %s variant.' % (sequence.id, sequence.lineage, cur_variant) )
                filtered_sequences.append(index)        
            index += 1
        return [sequences[i] for i in filtered_sequences]
                
    filtered_sequences = [] #stores the sequences that should be retained
    feasible_amplicons = set() #stores which amplicons are feasible
    
    lb = 0
    ub = 10**9
    
    #Filter based on the number of sequences of a lineage
    if min_sequences_threshold > 0 and lineages_location:
        sequences = filter_on_occurrences(sequences, lineages_location, min_sequences_threshold)
        
    #Filter based on variants
    if len(variants) > 0 and variants_location:
        sequences = filter_on_variants(sequences, variants_location, variants)
    
    #Determine feasible amplicons and indices where sequences potentially differ
    index = 0 #<- sequence index
    options_table = [set(['a','c','g','t','-']) for _ in range(sequences[0].length)] #This will contain sets with the possible nucleotides at every location in the sequence

    for sequence in sequences:
        sequence.align_to_trim()
        (cur_lb, cur_ub) = sequence.find_bounds(min_non_align)
        lb = max(lb, cur_lb)
        ub = min(ub, cur_ub)

        #Check the characters:
        #   By taking intersections we end up with a set that is either empty, or contains some character(s).
        #   In the former case, we cannot conclude whether all sequences are equal at an index, but in the
        #   latter case we know that there is a nucleotide that is shared by all sequences, hence allowing
        #   for a more efficient determination of amplicon differentiation.
        for c in range(sequence.length):
            options_table[c] = options_table[c].intersection(equivalent_characters(sequence.sequence[c]))
        
        if amplicon_width > 0 and misalign_threshold >= 0:
            if index == 0:
                feasible_amplicons = determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold)
            else:
                feasible_amplicons = feasible_amplicons.intersection(determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold))
        index += 1

    #Generate array with ones if an index should be considered (i.e. it isn't equivalent for all sequences and is contained in some feasible amplicon)
    relevant_nucleotides = np.zeros((sequences[0].length), dtype=np.int32)
    for amplicon in feasible_amplicons:
        for index in range(amplicon[0], amplicon[1]):
            if len(options_table[index]) == 0:
                relevant_nucleotides[index] = 1
    relevant_nucleotides = np.where(relevant_nucleotides == 1)[0]
    
    return sequences, lb, ub, feasible_amplicons, relevant_nucleotides

def generate_amplicons_mp_hybrid(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None, processors=1):
    '''
    Function that determines which sequence pairs can be differentiated for every amplicon in either $feasible_amplicons or in all possible amplicons.

    Parameters
    ----------
    sequences : list[ Sequence ]
        List of sequences that will be differentiated.
    amplicon_width : int
        Width of amplicons in number of nucleotides.
    comparison_matrix : dict [ (char,char) ]
        Dictionary that determines which characters should be considered equal.
    lb : int, optional
        Index from which amplicons should be generated if no feasible amplicons are supplied. The default is None in which case it is set to 0.
    ub : int, optional
        Last index (exclusive) where the final amplicon ends if they need to be generated. The default is None it is set to the length of the sequences.
    amplicon_threshold : int, optional
        Maximum number of allowed nucleotide mismatches between sequences in an amplicon. The default is 1.
    feasible_amplicons : set( (int,int) ), optional
        Set of amplicons which are defined as (start_index, end_index) where the end index is exclusive. The default is set() in which case amplicons will be generated.
    relevant_nucleotides : np.array, optional
        Numpy array with indices of nucleotides that can be different among sequences. The default is None.
    processors : int, optional
        Number of processors to use for multiprocessing. The default is 1.

    Returns
    -------
    res : list[ Amplicon ]
        List of final amplicons.
    X : np.array
        Numpy array containing 3 axes:
            amplicon
            sequence
            sequence
        where X[k,i,j] = 1 iff sequence i and j can be differentiated according to amplicon k.

    '''
    if not lb:
        lb = 0
    else:
        lb = max(lb,0)
    if not ub:
        ub = sequences[0].length
    else:
        ub = min(ub, sequences[0].length)
    
    #Check if feasible amplicons are provided
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        amplicons.sort(key = lambda x : x[0])
    else:
        amplicons = np.arange(lb, ub-lb-amplicon_width+1, 1, dtype=np.int32)
    
    #Transform input to numeric and otherwise relevant variables
    st = time.time()
    lineages_list = [sequence.lineage_num for sequence in sequences]
    ids_list = [sequence.id_num for sequence in sequences]
    _, comparison_matrix_num, sequences_num, AMPS = translate_to_numeric(sequences, amplicons, relevant_nucleotides, comparison_matrix)

    sequence_pairs_list = []
    for s1 in range(len(sequences)):
        for s2 in range(s1):
            if sequences[s1].lineage != sequences[s2].lineage:
                sequence_pairs_list.append([s2,s1])

    sequence_pairs_list = np.array(sequence_pairs_list, dtype=np.int32)
    ids_list = np.array(ids_list, dtype=np.int32)
    sequence_pairs_partition = [ sequence_pairs_list[i:i+ceil(sequence_pairs_list.shape[0]/processors)][:] for i in range(0, sequence_pairs_list.shape[0], ceil(sequence_pairs_list.shape[0]/processors)) ]
    partition_sizes = [spp.shape[0] for spp in sequence_pairs_partition]

    shared_mem = shared_memory.SharedMemory(create=True, size=np.dtype(np.int8).itemsize*np.prod((AMPS.shape[0], len(sequences), len(sequences))))
    shared_array = np.ndarray((len(amplicons), len(sequences), len(sequences)), dtype=np.int8, buffer=shared_mem.buf)
    print(str(len(amplicons)) + ', ' + str(len(sequences)))
    shared_array[:] = 0

    print(str(time.time() - st) + 's spent preprocessing')

    st = time.time()
    with mp.Pool(processors) as pool:
        pool.starmap(AmpliconGeneration.generate_amplicons_hybrid_cy, zip(itertools.repeat(AMPS), itertools.repeat(amplicon_width),
                                                                           itertools.repeat(AMPS.shape[0]), itertools.repeat(sequences_num),
                                                                           sequence_pairs_partition, partition_sizes, itertools.repeat(sequences_num.shape[0]),
                                                                           itertools.repeat(ids_list), itertools.repeat(comparison_matrix_num),
                                                                           itertools.repeat(relevant_nucleotides), itertools.repeat(relevant_nucleotides.shape[0]),
                                                                           itertools.repeat(amplicon_threshold), itertools.repeat(shared_mem.name)))
    print(str(time.time() - st) + 's spent differentiating')

    st = time.time()
    res = []
    for amplicon_index in range(AMPS.shape[0]):
        res.append(Amplicon(AMPS[amplicon_index][0], AMPS[amplicon_index][0]+amplicon_width))
    print(str(time.time() - st) + 's spent combining')

    st = time.time()
    X = np.zeros(shared_array.shape, dtype=np.int8)
    X[:] = shared_array[:]
    print(str(time.time() - st) + 's spent copying array')

    shared_mem.close()
    shared_mem.unlink()

    return res, X

def generate_amplicons_sp_hybrid(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None):
    '''
    Function that determines which sequence pairs can be differentiated for every amplicon in either $feasible_amplicons or in all possible amplicons.

    Parameters
    ----------
    sequences : list[ Sequence ]
        List of sequences that will be differentiated.
    amplicon_width : int
        Width of amplicons in number of nucleotides.
    comparison_matrix : dict [ (char,char) ]
        Dictionary that determines which characters should be considered equal.
    lb : int, optional
        Index from which amplicons should be generated if no feasible amplicons are supplied. The default is None in which case it is set to 0.
    ub : int, optional
        Last index (exclusive) where the final amplicon ends if they need to be generated. The default is None it is set to the length of the sequences.
    amplicon_threshold : int, optional
        Maximum number of allowed nucleotide mismatches between sequences in an amplicon. The default is 1.
    feasible_amplicons : set( (int,int) ), optional
        Set of amplicons which are defined as (start_index, end_index) where the end index is exclusive. The default is set() in which case amplicons will be generated.
    relevant_nucleotides : np.array, optional
        Numpy array with indices of nucleotides that can be different among sequences. The default is None.
    processors : int, optional
        Number of processors to use for multiprocessing. The default is 1.

    Returns
    -------
    res : list[ Amplicon ]
        List of final amplicons.
    X : np.array
        Numpy array containing 3 axes:
            amplicon
            sequence
            sequence
        where X[k,i,j] = 1 iff sequence i and j can be differentiated according to amplicon k.

    '''
    if not lb:
        lb = 0
    else:
        lb = max(lb,0)
    if not ub:
        ub = sequences[0].length
    else:
        ub = min(ub, sequences[0].length)
    
    #Check if feasible amplicons are provided
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        amplicons.sort(key = lambda x : x[0])
    else:
        amplicons = np.arange(lb, ub-lb-amplicon_width+1, 1, dtype=np.int32)
    
    #Transform input to numeric and otherwise relevant variables
    st = time.time()
    lineages_list = [sequence.lineage_num for sequence in sequences]
    ids_list = [sequence.id_num for sequence in sequences]
    _, comparison_matrix_num, sequences_num, AMPS = translate_to_numeric(sequences, amplicons, relevant_nucleotides, comparison_matrix)

    sequence_pairs_list = []
    for s1 in range(len(sequences)):
        for s2 in range(s1):
            if sequences[s1].lineage != sequences[s2].lineage:
                sequence_pairs_list.append([s2,s1])

    sequence_pairs_list = np.array(sequence_pairs_list, dtype=np.int32)
    ids_list = np.array(ids_list, dtype=np.int32)

    print(str(len(amplicons)) + ', ' + str(len(sequences)))
    print(str(time.time() - st) + 's spent preprocessing')

    st = time.time()
    X = AmpliconGeneration.generate_amplicons_hybrid_sp_cy(AMPS, amplicon_width, AMPS.shape[0], sequences_num, sequence_pairs_list, len(sequence_pairs_list), sequences_num.shape[0],
                                                            ids_list, comparison_matrix_num, relevant_nucleotides, relevant_nucleotides.shape[0], amplicon_threshold)
    X = np.asarray(X, dtype=np.int8)
    print(str(time.time() - st) + 's spent differentiating')

    st = time.time()
    res = []
    for amplicon_index in range(AMPS.shape[0]):
        res.append(Amplicon(AMPS[amplicon_index][0], AMPS[amplicon_index][0]+amplicon_width))
    print(str(time.time() - st) + 's spent combining')

    return res, X

def calculate_differences_per_amplicon_mp(sequences, amplicon_width, comparison_matrix, amplicon_threshold=1, processors=1):
    '''
    Function that determines which sequence pairs can be differentiated for every amplicon in either $feasible_amplicons or in all possible amplicons.

    Parameters
    ----------
    sequences : list[ Sequence ]
        List of sequences that will be differentiated.
    amplicon_width : int
        Width of amplicons in number of nucleotides.
    comparison_matrix : dict [ (char,char) ]
        Dictionary that determines which characters should be considered equal.
    amplicon_threshold : int, optional
        Maximum number of allowed nucleotide mismatches between sequences in an amplicon. The default is 1.
    processors : int, optional
        Number of processors to use for multiprocessing. The default is 1.

    Returns
    -------
    res : np.array
        Numpy array with the differences per amplicon (start index is equal to index in array).

    '''
    #Transform input to numeric and otherwise relevant variables
    ids_list = [sequence.id_num for sequence in sequences]
    _, comparison_matrix_num, sequences_num, _ = translate_to_numeric(sequences, [], [], comparison_matrix)

    sequence_pairs_list = []
    for s1 in range(len(sequences)):
        for s2 in range(s1):
            if sequences[s1].lineage != sequences[s2].lineage:
                sequence_pairs_list.append([s2,s1])
    max_possible_differences = len(sequence_pairs_list)

    sequence_pairs_list = np.array(sequence_pairs_list, dtype=np.int32)
    sequence_pairs_partition = [ sequence_pairs_list[i:i+ceil(sequence_pairs_list.shape[0]/processors)][:] for i in range(0, sequence_pairs_list.shape[0], ceil(sequence_pairs_list.shape[0]/processors)) ]
    partition_sizes = [spp.shape[0] for spp in sequence_pairs_partition]
    num_amplicons = sequences[0].length - amplicon_width + 1
    print(num_amplicons)

    with mp.Pool(processors) as pool:
        diffs = pool.starmap(AmpliconGeneration.calculate_amplicon_differences_cy, zip(itertools.repeat(amplicon_width), itertools.repeat(num_amplicons), 
                                                                            itertools.repeat(sequences_num), itertools.repeat(sequences[0].length),
                                                                            sequence_pairs_partition, partition_sizes, itertools.repeat(len(sequences)),
                                                                            itertools.repeat(comparison_matrix_num), itertools.repeat(amplicon_threshold)))
    res = diffs[0]
    for i in range(1, len(diffs)):
        res += diffs[i]
    return res, max_possible_differences

def translate_to_numeric(sequences, amplicons, relevant_nucleotides, comparison_matrix):
    
    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']
    char_comp = np.zeros((len(chars), len(chars)), dtype=np.int8)
    chars2num = {}
    seqs_num = np.zeros((len(sequences), sequences[0].length), dtype=np.int8)
    AMPS = np.zeros((len(amplicons), 3), dtype=np.int32)
    
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
        AMPS[a][0] = amplicons[a][0]
        cur = np.where(relevant_nucleotides < amplicons[a][0])[0]
        if cur.shape[0] > 0:
            AMPS[a][1] = cur[-1]
        cur = np.where(relevant_nucleotides < amplicons[a][1])[0]
        if cur.shape[0] > 0:
            AMPS[a][2] = cur[-1]
                
    return chars2num, char_comp, seqs_num, AMPS

def greedy(sequences, amplicons, differences_per_amplicon, primer_width, search_width, primer_index, comparison_matrix, max_amplicons, coverage, temperature_range, logging=False, multiplex=False):
    to_cover = np.sum(differences_per_amplicon, axis=0)
    to_cover = np.sum(to_cover > 0)

    if logging:
        log_results = ['Total to cover based on amplicon feasibility: ' + str(to_cover) + ' with ' + str(len(sequences)) + ' sequences and ' + str(len(amplicons)) + ' amplicons.']
    
    result_amplicons = []       #this will store the result amplicons

    amplicons.sort(key = lambda x : np.sum(differences_per_amplicon[x.id_num]), reverse=True) #sort based on differentiability
    while to_cover > 0 and len(result_amplicons) < max_amplicons and len(amplicons) > 0:
        best_amplicon = amplicons.pop(0)

        if logging:
            log_results.append('Checking amplicon: ' + str(best_amplicon.id))
        primer_index.check_amplicon(sequences, best_amplicon, primer_width, search_width)
        result_amplicons.append(best_amplicon)

        #Check if current amplicon can be added based on primer feasibility
        if multiplex:
            check = check_primer_feasibility(sequences, result_amplicons, primer_index, optimize=0, temperature_range=temperature_range, coverage=coverage)
        else:
            check = check_primer_feasibility(sequences, [best_amplicon], primer_index, optimize=0, temperature_range=temperature_range, coverage=coverage)
        
        #If amplicon can be added
        if check:
            to_cover = to_cover - np.sum(differences_per_amplicon[best_amplicon.id_num])
            if logging:
                log_results.append('Amplicon ' + str(best_amplicon.id) + ' succesfully added, new sequence pairs covered: ' + str(np.sum(differences_per_amplicon[best_amplicon.id_num])))
            for amplicon in amplicons:
                differences_per_amplicon[amplicon.id_num][differences_per_amplicon[best_amplicon.id_num] == 1] = 0
            #amplicons = [a for a in amplicons if np.sum(differences_per_amplicon[a.id_num]) > 0 and abs(a.start - best_amplicon.start) >= (a.end - a.start) + search_width]
            amplicons = [a for a in amplicons if np.sum(differences_per_amplicon[a.id_num]) > 0]
            amplicons.sort(key = lambda x : np.sum(differences_per_amplicon[x.id_num]), reverse=True)
        else:
            result_amplicons.pop(-1)
            if logging:
                log_results.append('Amplicon ' + str(best_amplicon.id) + ' rejected due to being unable to find primers to cover enough sequences')
    if logging:
        return log_results, result_amplicons
    else:
        return result_amplicons

def greedy_fancy(sequences, amplicons, differences_per_amplicon, primer_width, search_width, primer_index, comparison_matrix, max_amplicons, coverage, temperature_range, logging=False, multiplex=False):
    to_cover = np.sum(differences_per_amplicon, axis=0)
    to_cover = np.sum(to_cover > 0)

    if logging:
        log_results = ['Total to cover based on amplicon feasibility: ' + str(to_cover) + ' with ' + str(len(sequences)) + ' sequences and ' + str(len(amplicons)) + ' amplicons.']
    
    result_amplicons = []       #this will store the result amplicons
    result_primers = []         #this will store the result primers

    amplicons.sort(key = lambda x : np.sum(differences_per_amplicon[x.id_num]), reverse=True) #sort based on differentiability
    while to_cover > 0 and len(result_amplicons) < max_amplicons and len(amplicons) > 0:
        best_amplicon = amplicons.pop(0)

        if logging:
            log_results.append('Checking amplicon: ' + str(best_amplicon.id))
        primer_index.check_amplicon(sequences, best_amplicon, primer_width, search_width)
        result_amplicons.append(best_amplicon)

        #Check if current amplicon can be added based on primer feasibility
        [check, cur_primers, covered_differences, sequences_covered] = check_primer_feasibility_single_amplicon_max_coverage(sequences, best_amplicon, differences_per_amplicon[best_amplicon.id_num], np.sum(differences_per_amplicon[best_amplicon.id_num]), primer_index, temperature_range=temperature_range, coverage=coverage)
        if check:
            to_cover = to_cover - np.sum(covered_differences)
            if logging:
                log_results.append('Amplicon ' + str(best_amplicon.id) + ' succesfully added, new sequence pairs covered: ' + str(np.sum(covered_differences)) + '(fraction differences covered: ' + str(np.sum(covered_differences)/np.sum(differences_per_amplicon[best_amplicon.id_num])) + '), (fraction sequences covered: ' + str(sequences_covered) + ')')
            for amplicon in amplicons:
                differences_per_amplicon[amplicon.id_num][covered_differences == 1] = 0
            amplicons = [a for a in amplicons if np.sum(differences_per_amplicon[a.id_num]) > 0]
            amplicons.sort(key = lambda x : np.sum(differences_per_amplicon[x.id_num]), reverse=True)
            result_primers.append(cur_primers)
        else:
            result_amplicons.pop(-1)
            if logging:
                log_results.append('Amplicon ' + str(best_amplicon.id) + ' rejected due to being unable to find primers to cover enough sequences')
    if logging:
        return log_results, result_amplicons, result_primers
    else:
        return result_amplicons

def check_primer_feasibility(sequences, amplicons, primer_index, optimize=0, temperature_range=5, coverage=1):
    env = gp.Env(empty=True)
    env.setParam('OutputFlag',0) #turns off logging
    env.start()
    
    model = gp.Model(env=env)
    model.ModelSense = GRB.MINIMIZE

    #Primer variables
    forward_primers = {} #primer -> (variable, temp)
    reverse_primers = {} #primer -> (variable, temp)
    covered_binary = {} #(sequence, amplicon) -> variable
    #Initialize variables
    for amplicon in amplicons:
        for sequence in amplicon.primers['forward']:
            for primer in amplicon.primers['forward'][sequence]:
                if primer not in forward_primers:
                    forward_primers[primer] = (model.addVar(vtype=GRB.BINARY, obj=optimize), primer_index.index2primer['forward'][primer].temperature)
            for primer in amplicon.primers['reverse'][sequence]:
                if primer not in reverse_primers:
                    reverse_primers[primer] = (model.addVar(vtype=GRB.BINARY, obj=optimize), primer_index.index2primer['reverse'][primer].temperature)
            covered_binary[(sequence, amplicon.id)] = model.addVar(vtype=GRB.BINARY, obj=0)
    #Temperature variables
    max_temp = model.addVar(vtype=GRB.CONTINUOUS, obj=0)
    min_temp = model.addVar(vtype=GRB.CONTINUOUS, obj=0)
            
    #Enforce covered_binary variables to 0 if they aren't covered
    for amplicon in amplicons:
        #Coverage per sequence per amplicon
        for sequence in sequences:
            model.addConstr(covered_binary[(sequence.id_num, amplicon.id)] <= sum(forward_primers[primer][0] for primer in amplicon.primers['forward'][sequence.id_num]))
            model.addConstr(covered_binary[(sequence.id_num, amplicon.id)] <= sum(reverse_primers[primer][0] for primer in amplicon.primers['reverse'][sequence.id_num]))
        #At least $coverage (fraction) of the sequences should be covered per amplicon
        model.addConstr(sum(covered_binary[(sequence.id_num, amplicon.id)] for sequence in sequences) >= coverage * len(sequences))
        #Temperature constraints
        for primer in amplicon.full_primerset['forward']: #iterate over forward primers
            model.addConstr( min_temp <= primer_index.index2primer['forward'][primer].temperature * (3 - 2 * forward_primers[primer][0]) )
        for primer in amplicon.full_primerset['reverse']:
            model.addConstr( max_temp >= primer_index.index2primer['reverse'][primer].temperature * reverse_primers[primer][0] )
    model.addConstr(max_temp - min_temp <= temperature_range)
        
    #Check primer conflicts
    for pair in itertools.combinations(forward_primers.keys(), 2):
        model.addConstr( forward_primers[pair[0]][0] + forward_primers[pair[1]][0] <= primer_index.check_conflict( [primer_index.index2primer['forward'][pair[0]], primer_index.index2primer['forward'][pair[1]]] ) )
    for pair in itertools.combinations(reverse_primers.keys(), 2):
        model.addConstr( reverse_primers[pair[0]][0] + reverse_primers[pair[1]][0] <= primer_index.check_conflict( [primer_index.index2primer['reverse'][pair[1]], primer_index.index2primer['reverse'][pair[1]]] ) )
    for fwd in forward_primers:
        for rev in reverse_primers:
            model.addConstr( forward_primers[fwd][0] + reverse_primers[rev][0] <= primer_index.check_conflict( [primer_index.index2primer['forward'][fwd], primer_index.index2primer['reverse'][rev]] ) )
    
    model.optimize()
    if optimize == 1 and model.Status == 2:
        res = {'forward' : [], 'reverse' : []}
        print('Forward primers: ')
        for primer in forward_primers:
            if forward_primers[primer][0].x > 0.9:
                print(primer_index.index2primer['forward'][primer].sequence)
                res['forward'].append(primer_index.index2primer['forward'][primer].sequence)
        print('Reverse primers: ')
        for primer in reverse_primers:
            if reverse_primers[primer][0].x > 0.9:
                print(primer_index.index2primer['reverse'][primer].sequence)
                res['reverse'].append(primer_index.index2primer['reverse'][primer].sequence)
        return res
    return model.Status == 2

def check_primer_feasibility_single_amplicon_max_coverage(sequences, amplicon, differences, total_differences, primer_index, temperature_range=5, coverage=1):
    env = gp.Env(empty=True)
    env.setParam('OutputFlag',0)
    env.start()

    model = gp.Model(env=env)
    model.ModelSense = GRB.MAXIMIZE

    #Primer variables
    forward_primers = {} #primer -> (variable, temperature)
    reverse_primers = {} #primer -> (variable, temperature)
    #Sequence variables
    covered_binary = {} #sequence_id -> variable
    covered_pairs = {} #(sequence_id, sequence_id) -> variable

    #Initialize primer and sequence variables
    for sequence in amplicon.primers['forward']:
        for primer in amplicon.primers['forward'][sequence]:
            if primer not in forward_primers:
                forward_primers[primer] = (model.addVar(vtype=GRB.BINARY, obj=0), primer_index.index2primer['forward'][primer].temperature)
        for primer in amplicon.primers['reverse'][sequence]:
            if primer not in reverse_primers:
                reverse_primers[primer] = (model.addVar(vtype=GRB.BINARY, obj=0), primer_index.index2primer['reverse'][primer].temperature)
        covered_binary[sequence] = model.addVar(vtype=GRB.BINARY, obj=0)
    for s1 in range(len(sequences)):
        for s2 in range(s1):
            if sequences[s1].lineage != sequences[s2].lineage and differences[sequences[s2].id_num, sequences[s1].id_num] == 1:
                covered_pairs[(sequences[s1].id_num, sequences[s2].id_num)] = model.addVar(vtype=GRB.BINARY, obj=1)
                model.addConstr(covered_pairs[(sequences[s1].id_num, sequences[s2].id_num)] <= 0.5*covered_binary[sequences[s1].id_num] + 0.5*covered_binary[sequences[s2].id_num])
    num_primer_pairs = model.addVar(vtype=GRB.INTEGER, obj=-0.1*total_differences)

    #Temperature variables
    max_temp = model.addVar(vtype=GRB.CONTINUOUS, obj=0)
    min_temp = model.addVar(vtype=GRB.CONTINUOUS, obj=0)

    #Enforce covered_binary variables to 0 if they aren't covered
    for sequence in sequences:
        model.addConstr(covered_binary[sequence.id_num] <= sum(forward_primers[primer][0] for primer in amplicon.primers['forward'][sequence.id_num]))
        model.addConstr(covered_binary[sequence.id_num] <= sum(reverse_primers[primer][0] for primer in amplicon.primers['reverse'][sequence.id_num]))
        #At least $coverage (fraction) of the sequences should be covered per amplicon
        model.addConstr(sum(covered_binary[sequence.id_num] for sequence in sequences) >= coverage * len(sequences))
        #Temperature constraints
        for primer in amplicon.full_primerset['forward']: #iterate over forward primers
            model.addConstr( min_temp <= primer_index.index2primer['forward'][primer].temperature * (3 - 2 * forward_primers[primer][0]) )
        for primer in amplicon.full_primerset['reverse']:
            model.addConstr( max_temp >= primer_index.index2primer['reverse'][primer].temperature * reverse_primers[primer][0] )
    model.addConstr(max_temp - min_temp <= temperature_range)

    #Check primer conflicts
    for pair in itertools.combinations(forward_primers.keys(), 2):
        model.addConstr( forward_primers[pair[0]][0] + forward_primers[pair[1]][0] <= primer_index.check_conflict( [primer_index.index2primer['forward'][pair[0]], primer_index.index2primer['forward'][pair[1]]] ) )
    for pair in itertools.combinations(reverse_primers.keys(), 2):
        model.addConstr( reverse_primers[pair[0]][0] + reverse_primers[pair[1]][0] <= primer_index.check_conflict( [primer_index.index2primer['reverse'][pair[1]], primer_index.index2primer['reverse'][pair[1]]] ) )
    for fwd in forward_primers:
        for rev in reverse_primers:
            model.addConstr( forward_primers[fwd][0] + reverse_primers[rev][0] <= primer_index.check_conflict( [primer_index.index2primer['forward'][fwd], primer_index.index2primer['reverse'][rev]] ) )

    #Set variable for primer pairs
    model.addConstr(num_primer_pairs >= sum(forward_primers[primer][0] for primer in forward_primers))
    model.addConstr(num_primer_pairs >= sum(reverse_primers[primer][0] for primer in reverse_primers))

    model.optimize()
    if model.Status == 2:
        res = {'forward' : [], 'reverse' : []}
        seqs_covered = 0
        print('Forward primers: ')
        for primer in forward_primers:
            if forward_primers[primer][0].x > 0.9:
                print(primer_index.index2primer['forward'][primer].sequence)
                res['forward'].append(primer_index.index2primer['forward'][primer].sequence)
        print('Reverse primers: ')
        for primer in reverse_primers:
            if reverse_primers[primer][0].x > 0.9:
                print(primer_index.index2primer['reverse'][primer].sequence)
                res['reverse'].append(primer_index.index2primer['reverse'][primer].sequence)
        realized_differences = np.zeros(differences.shape, dtype=np.int8)
        for pair in covered_pairs:
            if covered_pairs[pair].x > 0.9:
                realized_differences[pair[1], pair[0]] = 1
        for sequence in covered_binary:
            if covered_binary[sequence].x > 0.9:
                seqs_covered += 1/len(sequences)
        return [True, res, realized_differences, seqs_covered]
    return [False, None, None, None]