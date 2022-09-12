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
        res[-1].differences = set(zip(cur_diffs[0], cur_diffs[1]))

    shared_mem.close()
    shared_mem.unlink()

    return res

#Smart and new and better
def generate_amplicons_mp_smart(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None, processors=1):
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
    #Dictionary that will contain differences per amplicon
    manager = mp.Manager()
    diffs_per_amplicon = manager.dict()
    for a in amplicons:
        diffs_per_amplicon[a] = '{'

    #Store information of interest in lists as to not pass entire sequence objects to Pool
    lineages_list = [sequence.lineage_num for sequence in sequences]
    sequences_list = [sequence.sequence for sequence in sequences]
    ids_list = [sequence.id_num for sequence in sequences]
    index_list = list(range(len(sequences)))

    sequence_pairs_list = list(itertools.combinations(index_list, 2))
    sequence_pair_partition = [ sequence_pairs_list[i:i+(ceil(len(sequence_pairs_list)/processors))] for i in range(0, len(sequence_pairs_list), ceil(len(sequence_pairs_list)/processors)) ]
    #amplicons_part = [ amplicons[i:i+(ceil(len(amplicons)/processors))] for i in range(0, len(amplicons), ceil(len(amplicons)/processors)) ]
    #included_amplicons = [ index_list_amplicons[i:i+(ceil(len(index_list_amplicons)/processors))] for i in range(0, len(index_list_amplicons), ceil(len(index_list_amplicons)/processors)) ]

    #Check which nucleotides have to be considered
    if not type(relevant_nucleotides) == np.ndarray:
        relevant_nucleotides = np.arange(lb,ub,1,dtype=int)

    with mp.Pool(processors) as pool:
        pool.starmap(AmpliconGeneration.determine_differences_smart_cy, 
        zip(itertools.repeat(amplicons), sequence_pair_partition, itertools.repeat(sequences_list), itertools.repeat(lineages_list), itertools.repeat(ids_list), 
        itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix), itertools.repeat(relevant_nucleotides),
        itertools.repeat(diffs_per_amplicon)))

    res = []
    for amplicon in amplicons:
        res.append(Amplicon(amplicon[0], amplicon[1]))
        res[-1].differences = eval(diffs_per_amplicon[amplicon][:len(diffs_per_amplicon[amplicon])-1] + '}')

    return res

def generate_amplicons_mp_smartest(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None, processors=1):
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
    #sequence_pairs_partition = [ sequence_pairs_list[i:i+ceil(len(sequence_pairs_list)/processors)] for i in range(0, len(sequence_pairs_list), ceil(len(sequence_pairs_list)/processors)) ]
    sequence_pairs_partition = [ sequence_pairs_list[i:i+ceil(sequence_pairs_list.shape[0]/processors)][:] for i in range(0, sequence_pairs_list.shape[0], ceil(sequence_pairs_list.shape[0]/processors)) ]
    partition_sizes = [spp.shape[0] for spp in sequence_pairs_partition]

    print(str(time.time() - st) + 's spent preprocessing')

    st = time.time()
    with mp.Pool(processors) as pool:
        """
        diffs_per_amp = pool.starmap(AmpliconGeneration.generate_amplicons_smarter_cy, zip(itertools.repeat(amplicons), itertools.repeat(amplicons_lb),
                                                                           itertools.repeat(amplicons_ub), itertools.repeat(amplicon_width),
                                                                           itertools.repeat(amplicons.shape[0]), itertools.repeat(sequences_num),
                                                                           sequence_pairs_partition, partition_sizes,
                                                                           itertools.repeat(ids_list), itertools.repeat(comparison_matrix_num),
                                                                           itertools.repeat(relevant_nucleotides), itertools.repeat(relevant_nucleotides.shape[0]),
                                                                           itertools.repeat(amplicon_threshold)))
        """
        diffs_per_amp = pool.starmap(AmpliconGeneration.generate_amplicons_smarter_cy, zip(itertools.repeat(AMPS), itertools.repeat(amplicon_width),
                                                                           itertools.repeat(AMPS.shape[0]), itertools.repeat(sequences_num),
                                                                           sequence_pairs_partition, partition_sizes,
                                                                           itertools.repeat(ids_list), itertools.repeat(comparison_matrix_num),
                                                                           itertools.repeat(relevant_nucleotides), itertools.repeat(relevant_nucleotides.shape[0]),
                                                                           itertools.repeat(amplicon_threshold)))
    print(str(time.time() - st) + 's spent differentiating')

    st = time.time()
    res = []
    for amplicon_index in range(AMPS.shape[0]):
        res.append(Amplicon(AMPS[amplicon_index][0], AMPS[amplicon_index][0]+amplicon_width))
        cur_diffs = [part[(AMPS[amplicon_index][0], AMPS[amplicon_index][0]+amplicon_width)] for part in diffs_per_amp]
        #res[-1].differences = set.union(*cur_diffs)
        res[-1].differences = set([item for sublist in cur_diffs for item in sublist])
    print(str(time.time() - st) + 's spent combining')
    print(len(res[0].differences))
    return res

def generate_amplicons_mp_hybrid(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), relevant_nucleotides=None, processors=1):
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
        """
        diffs_per_amp = pool.starmap(AmpliconGeneration.generate_amplicons_smarter_cy, zip(itertools.repeat(amplicons), itertools.repeat(amplicons_lb),
                                                                           itertools.repeat(amplicons_ub), itertools.repeat(amplicon_width),
                                                                           itertools.repeat(amplicons.shape[0]), itertools.repeat(sequences_num),
                                                                           sequence_pairs_partition, partition_sizes,
                                                                           itertools.repeat(ids_list), itertools.repeat(comparison_matrix_num),
                                                                           itertools.repeat(relevant_nucleotides), itertools.repeat(relevant_nucleotides.shape[0]),
                                                                           itertools.repeat(amplicon_threshold)))
        """
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
            model.addConstr(covered_binary[(sequence.alt_id, amplicon.id)] <= sum(forward_primers[primer][0] for primer in amplicon.primers['forward'][sequence.alt_id]))
            model.addConstr(covered_binary[(sequence.alt_id, amplicon.id)] <= sum(reverse_primers[primer][0] for primer in amplicon.primers['reverse'][sequence.alt_id]))
        #At least $coverage (fraction) of the sequences should be covered per amplicon
        model.addConstr(sum(covered_binary[(sequence.alt_id, amplicon.id)] for sequence in sequences) >= coverage * len(sequences))
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