from Bio import AlignIO
import csv
from math import ceil
import multiprocessing as mp
from multiprocessing import shared_memory
import itertools

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
                    sequences.append(Sequence(sequences_temp[meta[0]], meta[0], lineage=meta[skip]))
                
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

def generate_amplicons_mp_exp_cy(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), processors=1):

    #Check if there are feasible amplicons in input
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        amplicons.sort(key = lambda x : x[0])
    else:
        if not lb:
            lb = 0
        else:
            lb = max(lb, 0)
        if not ub:
            ub = sequences[0].length
        else:
            ub = min(ub, sequences[0].length)
        amplicons = [(i,i+amplicon_width) for i in range(ub - lb - amplicon_width + 1)]
    #Store information of interest in lists as to not pass entire sequence objects to Pool
    lins = [seq.lineage_num for seq in sequences]
    seqs = [seq.sequence for seq in sequences]
    ids = [seq.id_num for seq in sequences]
    #Partition the amplicon set which will get passed on different poolworkers
    amplicons_part = [ amplicons[i:i+(ceil(len(amplicons)/processors))] for i in range(0, len(amplicons), ceil(len(amplicons)/processors))]
    
    with mp.Pool(processors) as pool:
        res = pool.starmap(AmpliconGeneration.determine_differences_cy, zip(amplicons_part, itertools.repeat(seqs), itertools.repeat(lins), itertools.repeat(ids), itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix)))
    
    result_amplicons = []
    for partition in res:
        for amp in partition:
            result_amplicons.append(Amplicon(amp[0], amp[1]))
            result_amplicons[-1].differences = partition[amp]
    return result_amplicons    

def generate_amplicons_mp_exp_cy2(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None, processors=1):

    if not lb:
        lb = 0
    else:
        lb = max(lb, 0)
    if not ub:
        ub = sequences[0].length
    else:
        ub = min(ub, sequences[0].length)

    #Check if there are feasible amplicons in input
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        amplicons.sort(key = lambda x : x[0])
    else:
        amplicons = [(i,i+amplicon_width) for i in range(ub - lb - amplicon_width + 1)]
    #Store information of interest in lists as to not pass entire sequence objects to Pool
    lins = [seq.lineage_num for seq in sequences]
    seqs = [seq.sequence for seq in sequences]
    ids = [seq.id_num for seq in sequences]

    if not type(relevant_nucleotides) == np.ndarray:
        relevant_nucleotides = np.arange(lb, ub, 1, dtype = int)

    #Partition the amplicon set which will get passed on different poolworkers
    amplicons_part = [ amplicons[i:i+(ceil(len(amplicons)/processors))] for i in range(0, len(amplicons), ceil(len(amplicons)/processors))]
    to_check_part = [relevant_nucleotides[(relevant_nucleotides >= part[0][0])*(relevant_nucleotides <= part[-1][1])] for part in amplicons_part]
    
    with mp.Pool(processors) as pool:
        res = pool.starmap(AmpliconGeneration.determine_differences_cy2, zip(amplicons_part, itertools.repeat(seqs), itertools.repeat(lins), itertools.repeat(ids), itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix), to_check_part))

    result_amplicons = []
    for partition in res:
        for amp in partition:
            result_amplicons.append(Amplicon(amp[0], amp[1]))
            result_amplicons[-1].differences = partition[amp]
    return result_amplicons  

def generate_amplicons_mp_sequences(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None, processors=1):
    if not lb:
        lb = 0
    else:
        lb = max(lb, 0)
    if not ub:
        ub = sequences[0].length
    else:
        ub = min(ub, sequences[0].length)

    #Check for feasible amplicons in input
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        amplicons.sort(key = lambda x : x[0])
    else:
        amplicons = [(i,i+amplicon_width) for i in range(ub - lb - amplicon_width + 1)]

    #Store information of interest in lists as to not pass entire sequence objects to Pool
    lineages_list = [sequence.lineage_num for sequence in sequences]
    sequences_list = [sequence.sequence for sequence in sequences]
    ids_list = [sequence.id_num for sequence in sequences]
    index_list = list(range(len(sequences)))
    index_list_amplicons = list(range(len(amplicons)))

    sequence_pairs_list = list(itertools.combinations(index_list, 2))
    sequence_pair_partition = [ sequence_pairs_list[i:i+(ceil(len(sequence_pairs_list)/processors))] for i in range(0, len(sequence_pairs_list), ceil(len(sequence_pairs_list)/processors)) ]
    #amplicons_part = [ amplicons[i:i+(ceil(len(amplicons)/processors))] for i in range(0, len(amplicons), ceil(len(amplicons)/processors)) ]
    #included_amplicons = [ index_list_amplicons[i:i+(ceil(len(index_list_amplicons)/processors))] for i in range(0, len(index_list_amplicons), ceil(len(index_list_amplicons)/processors)) ]

    #Check which nucleotides have to be considered
    if not type(relevant_nucleotides) == np.ndarray:
        relevant_nucleotides = np.arange(lb,ub,1,dtype=int)

    shared_mem = shared_memory.SharedMemory(create=True, size=np.dtype(np.int8).itemsize*np.prod((len(amplicons), len(sequences), len(sequences))))
    shared_array = np.ndarray((len(amplicons), len(sequences), len(sequences)), dtype=np.int8, buffer=shared_mem.buf)
    shared_array[:] = 0

    with mp.Pool(processors) as pool:
        pool.starmap(AmpliconGeneration.determine_differences_sequencewise_cy, 
        zip(itertools.repeat(amplicons), sequence_pair_partition, itertools.repeat(sequences_list), itertools.repeat(lineages_list), itertools.repeat(ids_list), 
        itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix), itertools.repeat(relevant_nucleotides),
        itertools.repeat(shared_mem.name)))

    res = []
    for amplicon_index in range(len(amplicons)):
        res.append(Amplicon(amplicons[amplicon_index][0], amplicons[amplicon_index][1]))
        cur_diffs = np.where(shared_array[amplicon_index] == 1)
        res[-1].differences = res[-1].differences.union(set(zip(cur_diffs[0], cur_diffs[1])))

    shared_mem.close()
    shared_mem.unlink()

    return res