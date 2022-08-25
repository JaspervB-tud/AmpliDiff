from Bio import AlignIO
import RNA
from Bio.SeqUtils import MeltingTemp as mt
import os
import sys
import csv
import matplotlib.pyplot as plt
import itertools
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import copy
import time
import multiprocessing as mp
from Classes import *
#from pympler.asizeof import asizeof

#Updated scripts
'''############################################ This code deals with operations on str objects ##############################################################'''
def reverse_complement(sequence, rev=True):
    '''
    Function that returns the reverse complement of the given sequence

    Parameters
    ----------
    sequence : str
        String representation of the sequence to determine reverse complement of.
    rev : bool, optional
        True if the reverse should be returned, False otherwise. The default is True.

    Returns
    -------
    str
        Reverse complement of the input sequence.

    '''
    #Define the complement of every possible nucleotide
    translate = {
            'a' : 't',
            't' : 'a',
            'u' : 'a',
            'g' : 'c',
            'c' : 'g',
            'y' : 'r',
            'r' : 'y',
            's' : 's',
            'w' : 'w',
            'k' : 'm',
            'm' : 'k',
            'b' : 'v',
            'd' : 'h',
            'h' : 'd',
            'v' : 'b',
            'n' : 'n',
            '-' : '-'
        }
    res = ''
    for i in range(len(sequence)):
        res += translate[sequence[i]]
    if rev:
        return res[::-1]
    else:
        return res

def calculate_degeneracy(sequence):
    '''
    Function that returns the degeneracy of a sequence of nucleotides

    Parameters
    ----------
    sequence : str
        String representation of a series of consecutive nucleotides.

    Returns
    -------
    res : int
        Degeneracy of the input sequence.

    '''
    res = 1
    for char in sequence:
        if char in ['y', 'r', 's', 'w', 'm', 'k']:
            res = res*2
        elif char in ['b', 'd', 'h', 'v']:
            res = res*3
        elif char == 'n':
            res = res*4
    return res

def calculate_GC(sequence):
    '''
    Function that calculates the GC-content of a sequence

    Parameters
    ----------
    sequence : str
        String representation of sequence of interest.

    Returns
    -------
    float
        GC-content of input sequence.

    '''
    res = len(sequence)
    for char in sequence:
        if char in ['a','t','w','-']:
            res -= 1
    return res / len(sequence)

def calculate_end_stats(sequence, comparison_matrix):
    '''
    Function that determines the number of a/t (c/g) characters in the last 3 (5) characters of the input sequence

    Parameters
    ----------
    sequence : str
        String representation of sequence of interest.
    comparison_matrix : dict[ (char,char) ]
        Dictionary that determines which characters should be considered equal.

    Returns
    -------
    res : (int, int, bool)
        Triplet where the first element is the number of a/t chars in final 3, second element the number of c/g in final 5 and last element is true when last character is c/g.

    '''
    res = [0, 0, False]
    for i in range(1, 4):
        if comparison_matrix[(sequence[-i], 'a')][0] or comparison_matrix[(sequence[-i], 't')][0]:
            res[0] += 1
        elif comparison_matrix[(sequence[-i], 'c')][0] or comparison_matrix[(sequence[-i], 'g')][0]:
            res[1] += 1
            if i == 1:
                res[2] = True
    for i in range(4, 6):
        if comparison_matrix[(sequence[-i], 'c')][0] or comparison_matrix[(sequence[-i], 'g')][0]:
            res[1] += 1
    return res

def calculate_longest_monorun(sequence, comparison_matrix):
    '''
    Function that calculates the longest run of a single character in the given sequence

    Parameters
    ----------
    sequence : str
        String representation of sequence of interest.
    comparison_matrix : dict[ (char,char) ]
        Dictionary that determines which characters should be considered equal.

    Returns
    -------
    int
        Longest run of a single character.

    '''
    stats = [0, 1] # (longest, current)
    for i in range(1, len(sequence)):
        if comparison_matrix[(sequence[i-1], sequence[i])][0]:
            stats[1] += 1
        else:
            stats[0] = max(stats)
            stats[1] = 1
    return max(stats)

def calculate_longest_duorun(sequence, comparison_matrix):
    '''
    Function that calculates the longest run of a pair of characters in the given sequence

    Parameters
    ----------
    sequence : str
        String representation of sequence of interest.
    comparison_matrix : dict[ (char,char) ]
        Dictionary that determines which characters should be considered equal.

    Returns
    -------
    int
        Longest run of a pair of characters.

    '''
    stats = [0, 1] # (longest, current)
    current_duo = (sequence[0], sequence[1])
    index = 2
    while index < len(sequence) - 1:
        if comparison_matrix[(current_duo[0], sequence[index])][0] and comparison_matrix[(current_duo[1], sequence[index+1])][0]:
            stats[1] += 1
            index += 2
        else:
            stats[0] = max(stats)
            stats[1] = 1
            current_duo = (sequence[index-1], sequence[index])
            index += 1
    return max(stats)
            
def disambiguate(sequence):
    '''
    Function that disambiguates the given sequence by generating a string for every degenerate combination

    Parameters
    ----------
    sequence : str
        String representation of the sequence to disambiguate.

    Returns
    -------
    res : [str]
        List containing the non-degenerate sequences represented by the input sequence.

    '''
    translation = {'a' : ['a'],
                   'c' : ['c'],
                   'g' : ['g'],
                   't' : ['t'],                           
                   'b' : ['c','g','t'],
                   'd' : ['a','g','t'],
                   'h' : ['a','c','t'],
                   'k' : ['g','t'],
                   'm' : ['a','c'],
                   'n' : ['a','c','g','t'],
                   'r' : ['a','g'],
                   's' : ['g','c'],
                   'v' : ['a','c','g'],
                   'w' : ['a','t'],
                   'y' : ['c','t']}
    
    res = translation[sequence[0]].copy()
    for char in sequence[1:]:
        for subsequence_index in range(len(res)):
            new_subsequences = []
            for new_char in translation[char]:
                new_subsequences.append(res[subsequence_index] + new_char)
            res[subsequence_index] = new_subsequences
        res = list(itertools.chain(*res))
    return res
    
'''##########################################################################################################################################################'''

'''############################################ This code deals with generating generic objects #############################################################'''
def generate_opportunistic_matrix():
    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']
    char_comp = {
                'a' : ['a','r','m','w','d','h','v'],
                'c' : ['c','y','m','s','b','h','v'],
                't' : ['t','y','k','w','b','d','h'],
                'g' : ['g','r','k','s','b','d','v'],
                'u' : ['u','y','k','w','b','d','h'],
                'r' : ['a','g','r','k','m','s','w','b','d','h','v'],
                'y' : ['c','t','u','y','k','m','s','w','b','d','h','v'],
                'k' : ['g','t','u','r','y','k','s','w','b','d','h','v'],
                'm' : ['a','c','r','y','m','s','w','b','d','h','v'],
                's' : ['c','g','y','k','m','s','b','d','h','v'],
                'w' : ['a','t','u','r','y','k','m','w','b','d','h','v'],
                'b' : ['c','g','t','u','r','y','k','m','s','w','b','d','h','v'],
                'd' : ['a','g','t','u','r','y','k','m','s','w','b','d','h','v'],
                'h' : ['a','c','t','u','r','y','k','m','s','w','b','d','h','v'],
                'v' : ['a','c','g','r','y','k','m','s','w','b','d','h','v'],
                '-' : ['-']
        }
    res = {}
    for c1 in chars:
        for c2 in chars:
            if c1 == '-' or c2 == '-':
                if c1 == c2:
                    res[(c1,c2)] = (True,True)
                else:
                    res[(c1,c2)] = (False,True)
            elif c1 == 'n' or c2 == 'n':
                res[(c1,c2)] = (True,False)
            elif c1 in char_comp[c2] and c2 in char_comp[c1]:
                res[(c1,c2)] = (True,False)
            else:
                res[(c1,c2)] = (False,False)
    return res

'''##########################################################################################################################################################'''

'''############################################ This code deals with generating specific objects ############################################################'''
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
        if len(sequence.seq.replace('-','')) < 25000:
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

'''##########################################################################################################################################################'''

'''############################################ This code deals with general processing #####################################################################'''
def check_variant(sequence, lineages_per_variant, variants):
    '''
    Function that determines whether the sequence is part of one of the variants in $variants.

    Parameters
    ----------
    sequence : Sequence
        Sequence for which we want to know whether it is part of any of the variants of interest (not VOI).
    lineages_per_variant : dict[ str ]
        Dictionary that contains the lineages corresponding to every variant.
    variants : list[ str ]
        List with variants of interest (not VOI).

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
    $lb and $ub respectively. Additionally also determines feasible amplicons given amplicon width and misalignment character threshold.

    Parameters
    ----------
    sequences : list[ Sequence ]
        List with sequences to filter.
    min_non_align : int
        Number of nucleotides to include before (after) $lb ($ub).
    variants_location : str, optional
        Location of the Variants.txt file that explains which lineages are part of which variant. The default is None in which case every variant is considered.
    variants : list[ str ], optional
        List of the variants of interest (not VOI). The default is [] in which case every variant is considered.
    amplicon_width : int, optional
        Size of the amplicons, if we want to determine their feasibility a priori. The default is 0 in which case feasibility of amplicons is not checked.
    misalign_threshold : int, optional
        Number of allowed misalign characters in an amplicon. The default is -1 in which case feasibility of amplicons is not checked.
    lineages_location : str, optional
        Location of the lineages.txt file which contains the number of sequences per lineage. The default is None in which case every lineage is considered.
    min_sequences_threshold : int, optional
        Relative number of sequences that have to belong to a lineage in order to consider it. The default is 0 in which case every lineage is considered.

    Returns
    -------
    sequences : list[Sequence]
        List of sequences that belong the specified variants.
    lb : int
        Lowerbound such that every sequence has at least $min_non_align nucleotides before it.
    ub : int
        Upperbound such that every sequence has at least $min_non_align nucleotides after it.
    feasible_amplicons : set{ (int,int) }, optional
        Set of feasible amplicons in case amplicon feasibility should be checked.

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
                sequences_per_lineage[cur_line[0]] = int(cur_line[0])
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
                if sequences_per_lineage[sequence.lineage]/total_sequences >= min_sequences_threshold:
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
       
    index = 0 #<- sequence index
    for sequence in sequences:
        sequence.align_to_trim()
        (cur_lb, cur_ub) = sequence.find_bounds(min_non_align)
        lb = max(lb, cur_lb)
        ub = min(ub, cur_ub)
        
        if amplicon_width > 0 and misalign_threshold >= 0:
            if index == 0:
                feasible_amplicons = determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold)
            else:
                feasible_amplicons = feasible_amplicons.intersection(determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold))
        index += 1
    
    return sequences, lb, ub, feasible_amplicons
    
        
    """
    Old version
    
    #If we want to filter on the number of sequences
    if min_sequences_threshold > 0 and lineages_location:
        sequences_per_lineage = {}
        total_sequences = 0
        for line in csv.reader(open(lineages_location + '/lineages.txt'), delimiter='\t'):
            cur_line = line[0].split()
            #Tally the number of sequences of each lineage
            sequences_per_lineage[cur_line[0]] = int(cur_line[1])
            total_sequences += int(cur_line[1])
        print(total_sequences)
        #Check which sequences belong to a lineage that exceeds the minimum threshold
        index = 0
        for sequence in sequences:
            try:
                if sequences_per_lineage[sequence.lineage]/total_sequences >= min_sequences_threshold:
                    filtered_sequences.append(index)
            except:
                print(sequence.id + ' : ' + sequence.lineage + ' does not occur in lineages.txt')
            index += 1
        sequences = [sequences[i] for i in filtered_sequences]
        filtered_sequences = []
    
    if len(variants) > 0 and variants_location:
        lineages_per_variant = {} #Store the (base-)lineages per variant to infer variants later
        for variant in csv.reader(open(variants_location + '/Variants.txt'), delimiter='\t'):
            try:
                lineages_per_variant[variant[0]] = []
                for lineage in variant[1:]:
                    if lineage != '':
                        lineages_per_variant[variant[0]].append(lineage)
            except:
                continue
        index = 0
        first_sequence = True #set to false after finding the first feasible sequence
        for sequence in sequences:
            cur_variant = check_variant(sequence, lineages_per_variant, variants)
            if cur_variant:
                print( 'Sequence %s with lineage %s is part of the %s variant' % (sequence.id, sequence.lineage, cur_variant) )
                sequence.align_to_trim()
                (cur_lb, cur_ub) = sequence.find_bounds(min_non_align)
                lb = max(lb, cur_lb)
                ub = min(ub, cur_ub)
                filtered_sequences.append(index)
                #If amplicons should be checked, determine feasible amplicons for current sequence
                if amplicon_width > 0 and misalign_threshold >= 0:
                    if first_sequence:
                        feasible_amplicons = determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold)
                        first_sequence = False
                    else:
                        feasible_amplicons = feasible_amplicons.intersection(determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold))
            index += 1
        #If amplicons should be checked, return feasible amplicons
        if amplicon_width > 0 and misalign_threshold >= 0:
            return [sequences[i] for i in filtered_sequences], lb, ub, feasible_amplicons
        else:
            return [sequences[i] for i in filtered_sequences], lb, ub
    else:
        index = 0
        for sequence in sequences:
            sequence.align_to_trim()
            (cur_lb, cur_ub) = sequence.find_bounds(min_non_align)
            lb = max(lb, cur_lb)
            ub = min(ub, cur_ub)
            #If amplicons should be checked, determine feasible amplicons for current sequence
            if amplicon_width > 0 and misalign_threshold >= 0:
                if index == 0:
                    feasible_amplicons = determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold)
                else:
                    feasible_amplicons = feasible_amplicons.intersection(determine_feasible_amplicons(sequence, lb, ub, amplicon_width, misalign_threshold))
            index += 1
        #If amplicons should be checked, return feasible amplicons
        if amplicon_width > 0 and misalign_threshold >= 0:
            return sequences, lb, ub, feasible_amplicons
        else:
            return sequences, lb, ub
        
    """
def generate_amplicons(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=0, feasible_amplicons=set()):
    '''
    Function that generates amplicons, along with their differentiability. If a non-empty $feasible_amplicons set is provided, it will only
    generate the amplicons listed in $feasible_amplicons.
    For a multi-processing version of this see generate_amplicons_mp.

    Parameters
    ----------
    sequences : list[ Sequence ]
        List of sequences to base differentiability of amplicons on.
    amplicon_width : int
        Width of amplicons.
    comparison_matrix : dict[ (char,char) ]
        Dictionary that determines which characters should be considered equal.
    lb : int, optional
        Starting index  (inclusive) of amplicons. If feasible_amplicons is provided then this is ignored. The default is None in which case it is set to 0.
    ub : int, optional
        Ending index (exclusive) of amplicons. If feasible_amplicons is provided then this is ignored. The default is None in which case it is set to the sequence length.
    amplicon_threshold : int, optional
        Number of allowed mismatches between sequences. If two sequences differ at more than $amplicon_threshold indices then they are considered different. The default is 0.
    feasible_amplicons : set( (int,int) ), optional
        Set containing the amplicons (defined as start,end) for which differentiability has to be determined. The default is set() in which case all possible amplicons between lb and ub will be generated.

    Returns
    -------
    diffs_per_amp : dict [ (int,int) ] -> list [ (Sequence.id,Sequence.id) ]
        Dictionary containing the pairwise differences for every amplicon that has been considered.

    '''
    #Function that deals with moving the amplicon window to the next amplicon
    def move_window(seq1, seq2, diffs, prev_amplicon, next_amplicon):
        '''
        Auxilliary function that moves from the $prev_amplicon to the $next_amplicon by considering whether there is overlap or not.

        Parameters
        ----------
        seq1 : int
            Index of the first sequence to compare.
        seq2 : int
            Index of the second sequence to compare.
        diffs : list[ int ]
            List with indices (in previous amplicon) where sequences $seq1 and $seq2 differ.
        prev_amplicon : (int, int)
            Tuple with start (inclusive) and end (exclusive) indices of the previous amplicon.
        next_amplicon : (int, int)
            Tuple with start (inclusive) and end (exclusive) indices of the next amplicon.

        Returns
        -------
        diffs : list[ int ]
            List with indices where sequencces $seq1 and $seq2 differ in $next_amplicon

        '''
        #First check if there is overlap between the previous and next amplicon
        if prev_amplicon[1] > next_amplicon[0]: #overlap -> example: [2,12) [10,20) overlap of 2
            #Remove all indices that are not part of the next amplicon
            if len(diffs) > 0:
                while diffs[0] < next_amplicon[0]:
                    diffs.pop(0)
                    if len(diffs) == 0: #if the diffs list is empty, stop checking
                        break
            #Starting from the final index of the previous amplicon, find differences
            for i in range(prev_amplicon[1], next_amplicon[1]):
                if not comparison_matrix[(sequences[seq1].sequence[i], sequences[seq2].sequence[i])][0]: #check if next nucleotide matches
                    diffs.append(i)
        else: #no overlap
            diffs = []
            #Iterate over the next amplicon
            for i in range(next_amplicon[0], next_amplicon[1]):
                if not comparison_matrix[(sequences[seq1].sequence[i], sequences[seq2].sequence[i])][0]:
                    diffs.append(i)
        return diffs
    
    if len(feasible_amplicons) > 0: #set of feasible amplicons is provided
        diffs_per_amp = {a : None for a in feasible_amplicons}
    else: #check every amplicon
        #First check the lower and upperbounds
        if not lb:
            lb = 0
        else:
            lb = max(lb, 0)
        if not ub:
            ub = sequences[0].length
        else:
            ub = min(ub, sequences[0].length)
        diffs_per_amp = {(i, i+amplicon_width) : None for i in range(lb, ub-amplicon_width+1)}
    #Here diffs_per_amp contains all the necessary keys to generate the amplicons
    amplicons = list(diffs_per_amp.keys())
    amplicons.sort()
    
    for amplicon in amplicons:
        diffs_per_amp[amplicon] = [] #initialize empty list to save differentiation of amplicons
    for seq1 in range(len(sequences)):
        for seq2 in range(seq1):
            prev_amplicon = (0,0) #this makes sure that in the first call to move_window we actually determine the amplicon differentiability based on the first amplicon
            diffs = [] #start with an empty list of differences
            for next_amplicon in amplicons:
                diffs = move_window(seq1, seq2, diffs, prev_amplicon, next_amplicon) #determine differences in the next amplicon
                if len(diffs) > amplicon_threshold: #if the number of differences exceeds threshold, add sequence pair to diffs_per_amp
                    diffs_per_amp[next_amplicon].append((sequences[seq2].id, sequences[seq1].id))
                prev_amplicon = next_amplicon
    return diffs_per_amp 
            
'''############################################ This code deals with multiprocessing amplicon generation ####################################################'''
def move_window(seq1, seq2, diffs, prev_amplicon, next_amplicon, comparison_matrix):
    if prev_amplicon[1] > next_amplicon[0]:
        if len(diffs) > 0:
            while diffs[0] < next_amplicon[0]:
                diffs.pop(0)
                if len(diffs) == 0:
                    break
        for i in range(prev_amplicon[1], next_amplicon[1]):
            if not comparison_matrix[(seq1.sequence[i], seq2.sequence[i])][0]:
                diffs.append(i)
    else:
        diffs = []
        for i in range(next_amplicon[0], next_amplicon[1]):
            if not comparison_matrix[(seq1.sequence[i], seq2.sequence[i])][0]:
                diffs.append(i)
    return diffs

def determine_differences(amplicons, sequences, amplicon_threshold, comparison_matrix):
    diffs_per_amp = {a : set() for a in amplicons}
    
    for seq1 in range(len(sequences)):
        for seq2 in range(seq1):
            if sequences[seq1].lineage == sequences[seq2].lineage:
                pass
            else:
                prev_amplicon = (0,0)
                diffs = []
                for next_amplicon in amplicons:
                        diffs = move_window(sequences[seq1], sequences[seq2], diffs, prev_amplicon, next_amplicon, comparison_matrix)
                        if len(diffs) > amplicon_threshold:
                            diffs_per_amp[next_amplicon].add((sequences[seq2].id, sequences[seq1].id))
                        prev_amplicon = next_amplicon
    return diffs_per_amp
        
def generate_amplicons_mp2(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), processors=1, sort_type=0):
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        if sort_type == 0:
            amplicons.sort(key = lambda x : x[0])
        else:
            amplicons.sort()
    else:
        if not lb:
            lb = 0
        else:
            lb = max(lb, 0)
        if not ub:
            ub = sequences[0].length
        else:
            ub = min(ub, sequences[0].length)
        amplicons = [(i,i+amplicon_width) for i in range(lb, ub-amplicon_width+1)]
    amplicons_partitioned = [ amplicons[i:i+(ceil(len(amplicons)/processors))] for i in range(0, len(amplicons), ceil(len(amplicons)/processors))]

    with mp.Pool(processors) as pool:
        res = pool.starmap(determine_differences, zip(amplicons_partitioned, itertools.repeat(sequences), itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix)))
   
    result_amplicons = []
    for partition in res:
        for amp in partition:
            result_amplicons.append(Amplicon(amp[0], amp[1]))
            result_amplicons[-1].differences = partition[amp]
    return result_amplicons


''' Experimental!!
'''
def move_window_exp(seq1, seq2, diffs, prev_amplicon, next_amplicon, comparison_matrix):
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

def determine_differences_exp(amplicons, sequences, lineages, ids, amplicon_threshold, comparison_matrix):
    diffs_per_amp = {a : set() for a in amplicons}
    
    for seq1 in range(len(sequences)):
        for seq2 in range(seq1):
            if lineages[seq1] == lineages[seq2]:
                pass
            else:
                prev_amplicon = (0,0)
                diffs = []
                for next_amplicon in amplicons:
                    diffs = move_window_exp(sequences[seq1], sequences[seq2], diffs, prev_amplicon, next_amplicon, comparison_matrix)
                    if len(diffs) > amplicon_threshold:
                        diffs_per_amp[next_amplicon].add((ids[seq2], ids[seq1]))
                    prev_amplicon = next_amplicon
    return diffs_per_amp

def generate_amplicons_mp_exp(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), processors=1):
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
    lineages = [seq.lineage for seq in sequences]
    seqs = [seq.sequence for seq in sequences]
    ids = [seq.id for seq in sequences]
    amplicons_part = [ amplicons[i:i+(ceil(len(amplicons)/processors))] for i in range(0, len(amplicons), ceil(len(amplicons)/processors))]
    
    with mp.Pool(processors) as pool:
        res = pool.starmap(determine_differences_exp, zip(amplicons_part, itertools.repeat(seqs), itertools.repeat(lineages), itertools.repeat(ids), itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix)))
    
    result_amplicons = []
    for partition in res:
        for amp in partition:
            result_amplicons.append(Amplicon(amp[0], amp[1]))
            result_amplicons[-1].differences = partition[amp]
    return result_amplicons
'''
'''

def process_sequence_pairs(pairs, amplicons, amplicon_threshold, comparison_matrix):
    res = ((pair[0].id, pair[1].id), [])
    prev_amplicon = (0,0)
    diffs = []
    for amplicon in amplicons:
        diffs = move_window(pair[0], pair[1], diffs, prev_amplicon, amplicon, comparison_matrix)
        if len(diffs) > amplicon_threshold:
            res[1].append(amplicon)
        prev_amplicon = amplicon
    return res
                
def generate_amplicons_mp(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, feasible_amplicons=set(), processors=1):
    if len(feasible_amplicons) > 0:
        amplicons = list(feasible_amplicons)
        amplicons.sort()
    else:
        #First check the lower and upperbounds
        if not lb:
            lb = 0
        else:
            lb = max(lb, 0)
        if not ub:
            ub = sequences[0].length
        else:
            ub = min(ub, sequences[0].length)
            
        amplicons = [(i,i+amplicon_width) for i in range(lb, ub-amplicon_width+1)]
    #Filter out pairs of sequences of the same lineage
    sequence_pairs = []
    for pair in itertools.combinations(sequences, 2):
        if pair[0].lineage != pair[1].lineage:
            sequence_pairs.append(pair)
    
    sequence_pairs_partitioned = [ sequence_pairs[i:i+(ceil(len(sequence_pairs)/processors))] for i in range(0, len(sequence_pairs), ceil(len(sequence_pairs)/processors))]
    with mp.Pool(processors) as pool:
        check = pool.starmap(process_sequence_pair, zip(sequence_pairs, itertools.repeat(amplicons), itertools.repeat(amplicon_threshold), itertools.repeat(comparison_matrix)))
    return translate_amplicons(check)
    
def translate_amplicons(diffs_per_pair):
    amplicons = []
    for sequence_pair in diffs_per_pair:
        for amplicon in sequence_pair[1]:
            if amplicon not in amplicons:
                amplicons.append(Amplicon(amplicon[0], amplicon[1]))
                amplicons[-1].differences.add(sequence_pair[0])
            else:
                amplicons[amplicons.index(amplicon)].differences.add(sequence_pair[0])
    return amplicons
'''##########################################################################################################################################################'''

'''############################################ This code deals with the greedy amplicon selection ##########################################################'''
def greedy(sequences, amplicons, primer_width, search_width, primer_index, comparison_matrix, max_amplicons, coverage, temperature_range, logging=False, multiplex=True):
    to_cover = set()
    for amplicon in amplicons:
        to_cover = to_cover.union(amplicon.differences)
    total_to_cover = len(to_cover)
    
    if logging:
        log_results = ['Total to cover based on amplicon feasibility: ' + str(total_to_cover) + ' with ' + str(len(sequences)) + ' sequences and ' + str(len(amplicons)) + ' amplicons.']
        
    result_amplicons = [] #this will store the result amplicons
    removed_amplicons = [] #this will store the removed amplicons
    
    amplicons.sort(key = lambda x : len(x.differences), reverse=True) #sort the amplicons starting with the most differentiable amplicon
    while len(to_cover) > 0 and len(result_amplicons) < max_amplicons and len(amplicons) > 0:
        best_amplicon = amplicons.pop(0)
        removed_amplicons.append(best_amplicon)
        to_remove = [] #this will store the amplicons to remove
        
        if logging:
            log_results.append('Checking amplicon: ' + str(best_amplicon.id))   
        primer_index.check_amplicon(sequences, best_amplicon, primer_width, search_width)
        result_amplicons.append(best_amplicon)
        #Check if we can add current amplicon
        if multiplex:
            check = check_primer_feasibility(sequences, result_amplicons, primer_index, optimize=0, temperature_range=temperature_range, coverage=coverage)
        else:
            check = check_primer_feasibility(sequences, [best_amplicon], primer_index, optimize=0, temperature_range=temperature_range, coverage=coverage)
        #If amplicon can be added:    
        if check:
            to_cover = to_cover - best_amplicon.differences
            if logging:
                log_results.append('Amplicon ' + str(best_amplicon.id) + ' succesfully added, new sequence pairs covered: ' + str(len(best_amplicon.differences)))
            for amplicon in amplicons:
                amplicon.differences = amplicon.differences - best_amplicon.differences
                if len(amplicon.differences) == 0 or abs(amplicon.start - best_amplicon.start) < (amplicon.end - amplicon.start) + search_width:
                    to_remove.append(amplicon)
            for amplicon in to_remove:
                removed_amplicons.append(amplicon)
                amplicons.remove(amplicon)
            amplicons.sort(key = lambda x : len(x.differences), reverse=True)
        else:
            result_amplicons.pop(-1)
            if logging:
                log_results.append('Amplicon ' + str(best_amplicon.id) + ' rejected due to being unable to find primers to cover enough sequences')
    amplicons = amplicons + removed_amplicons
    if logging:
        return log_results, amplicons, result_amplicons
    else:
        return amplicons, result_amplicons
        
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
            model.addConstr(covered_binary[(sequence.id, amplicon.id)] <= sum(forward_primers[primer][0] for primer in amplicon.primers['forward'][sequence.id]))
            model.addConstr(covered_binary[(sequence.id, amplicon.id)] <= sum(reverse_primers[primer][0] for primer in amplicon.primers['reverse'][sequence.id]))
        #At least $coverage (fraction) of the sequences should be covered per amplicon
        model.addConstr(sum(covered_binary[(sequence.id, amplicon.id)] for sequence in sequences) >= coverage * len(sequences))
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


if __name__ == '__main__':
    sequences = generate_sequences('/Users/jaspervanbemmelen/Documents/Wastewater/Data/Global/global_all_time_N0','/Users/jaspervanbemmelen/Documents/Wastewater/Data/Global/global_all_time_N0')
    #sequences = generate_sequences('/Users/jaspervanbemmelen/Documents/Wastewater/Data/Africa/South_Africa/South_Africa_jan_2022', '/Users/jaspervanbemmelen/Documents/Wastewater/Data/Africa/South_Africa/South_Africa_jan_2022')
    sequences1, lb1, ub1, feasible1 = preprocess_sequences(sequences, 50, variants_location='/Users/jaspervanbemmelen/Documents/Wastewater/Data', amplicon_width=200, misalign_threshold=5)
    sequences2, lb2, ub2, feasible2 = preprocess_sequences(sequences, 50, variants_location='/Users/jaspervanbemmelen/Documents/Wastewater/Data', amplicon_width=200, misalign_threshold=5, lineages_location='/Users/jaspervanbemmelen/Documents/Wastewater/source_code', min_sequences_threshold=0.001)
    
    '''
    test = []
    lins = []
    for sequence in sequences:
        if sequence.lineage not in lins:
            test.append(sequence)
            lins.append(sequence.lineage)
    
    amplicon_width = 200
    comparison_matrix = generate_opportunistic_matrix()
    feasible_amplicons = feasible5
    amplicon_threshold = 1
    n_cores = 8
    n_seqs = 650
    
    #PI = PrimerIndex.generate_index_mp(sequences, 25, comparison_matrix, processors=3)
    
    st = time.time()
    amplicons2 = generate_amplicons_mp_exp(test[:n_seqs], amplicon_width, comparison_matrix, feasible_amplicons=feasible5, processors=n_cores)
    print(time.time() - st)
    
    st = time.time()
    amplicons = generate_amplicons_mp2(test[:n_seqs], amplicon_width, comparison_matrix, feasible_amplicons=feasible5, processors=n_cores, sort_type=0)
    print(time.time() - st)
    '''
"""
################################# This code generates the comparison matrices ######################################################
def generateConservativeMatrix():
    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']
    char_comp = {
                'a' : ['a','r','m','w','d','h','v'],
                'c' : ['c','y','m','s','b','h','v'],
                't' : ['t','y','k','w','b','d','h'],
                'g' : ['g','r','k','s','b','d','v'],
                'u' : ['u','y','k','w','b','d','h'],
                'r' : ['a','g','d','v'],
                'y' : ['c','t','u','b'],
                'k' : ['g','t','u','b','d'],
                'm' : ['a','c','h','v'],
                's' : ['c','g','b','h','v'],
                'w' : ['a','t','u'],
                'b' : ['c','g','t','u','y','k','s'],
                'd' : ['a','g','t','u','r','k','w'],
                'h' : ['a','c','t','u','y','m','w'],
                'v' : ['a','c','g','r','m','s'],
                '-' : ['-']
        }
    res = {}
    for c1 in chars:
        for c2 in chars:
            if c1 == '-' or c2 == '-':
                if c1 == c2:
                    res[(c1,c2)] = (True,True)
                else:
                    res[(c1,c2)] = (False,True)
            elif c1 == 'n' or c2 == 'n':
                res[(c1,c2)] = (True,False)
            elif c1 in char_comp[c2] and c2 in char_comp[c1]:
                res[(c1,c2)] = (True,False)
            else:
                res[(c1,c2)] = (False,False)
    return res

def generateOpportunisticMatrix():
    chars = ['a','c','t','g','u','r','y','k','m','s','w','b','d','h','v','n','-']
    char_comp = {
                'a' : ['a','r','m','w','d','h','v'],
                'c' : ['c','y','m','s','b','h','v'],
                't' : ['t','y','k','w','b','d','h'],
                'g' : ['g','r','k','s','b','d','v'],
                'u' : ['u','y','k','w','b','d','h'],
                'r' : ['a','g','r','k','m','s','w','b','d','h','v'],
                'y' : ['c','t','u','y','k','m','s','w','b','d','h','v'],
                'k' : ['g','t','u','r','y','k','s','w','b','d','h','v'],
                'm' : ['a','c','r','y','m','s','w','b','d','h','v'],
                's' : ['c','g','y','k','m','s','b','d','h','v'],
                'w' : ['a','t','u','r','y','k','m','w','b','d','h','v'],
                'b' : ['c','g','t','u','r','y','k','m','s','w','b','d','h','v'],
                'd' : ['a','g','t','u','r','y','k','m','s','w','b','d','h','v'],
                'h' : ['a','c','t','u','r','y','k','m','s','w','b','d','h','v'],
                'v' : ['a','c','g','r','y','k','m','s','w','b','d','h','v'],
                '-' : ['-']
        }
    res = {}
    for c1 in chars:
        for c2 in chars:
            if c1 == '-' or c2 == '-':
                if c1 == c2:
                    res[(c1,c2)] = (True,True)
                else:
                    res[(c1,c2)] = (False,True)
            elif c1 == 'n' or c2 == 'n':
                res[(c1,c2)] = (True,False)
            elif c1 in char_comp[c2] and c2 in char_comp[c1]:
                res[(c1,c2)] = (True,False)
            else:
                res[(c1,c2)] = (False,False)
    return res

################################# This code deals with preprocessing stuff ######################################################
def preprocessSequences(sequences, min_non_align):
    '''
    Function that preprocesses the input sequences to determine: a map from the aligned sequence to its raw form for every sequence, and
    the indices in the aligned sequence before which and after which there are at least $min_non_align non misalignment characters

    Parameters
    ----------
    sequences : list[Sequence]
        List of the sequences to be parsed.
    min_non_align : int
        Number of non-misalignment characters that should occur before lb, and after ub respectively.

    Returns
    -------
    lb : int
        Index before which there occur at least $min_non_align$ non-alignment characters.
    ub : int
        Index after which there occur at least $min_non_align$ non-alignment characters.

    '''
    lb = 0
    ub = sys.maxsize
    for sequence in sequences:
        sequence.alignToTrim()
        tmp1, tmp2 = sequence.findBounds(min_non_align)
        lb = max(lb, tmp1)
        ub = min(ub, tmp2)
        
    return lb, ub

def saveSequences(sequences, output_location, output_name):
    out_loc = os.path.abspath(output_location + output_name)
    with open(out_loc, 'w') as f:
        for sequence in sequences:
            f.write(sequence.toText() + '\n')

################################# This code generates and filters amplicons   #######################################################
def generateAmplicons(sequences, amplicon_width, comparison_matrix, lb=None, ub=None, amplicon_threshold=1, misalign_threshold=sys.maxsize):
    '''
    Function that generates a set of feasible amplicons for the given sequences.

    Parameters
    ----------
    sequences : list[Sequence]
        List of Sequence objects of equal length.
    amplicon_width : int
        Number of characters in the amplicons.
    comparison_matrix : dict[(char,char)]
        Comparison matrix which for every character pair returns a boolean tuple where the first element is true if chars are considered identical,
        and the second element is true if either character is a misalign character '-'
    lb : int, optional
        Lowerbound for the amplicon indices. The default is None.
    ub : int, optional
        Upperbound for the amplicon indices. The default is None.
    amplicon_threshold : int, optional
        Number of characters (nucleotides) that may differ at an amplicon in determining whether two sequences equivalent at an amplicon. The default is 1.
    misalign_threshold : int, optional
        Number of allowed misalign characters in amplicons. For example if s1 = act--ag and s2 = ac-ggag, then this amplicon is only viable
        if $misalign_threshold$ > 3. The default is Inf.

    Returns
    -------
    diffs_per_amp : dict[(int,int)]
        dictionary with amplicons (start,end(exclusive)) as keys and a dictionary with the sequences that differ for this amplicon.

    '''
    #Function that deals with moving the amplicon window 1 step to the right
    def moveWindow(D, mis, index, amplicon_width):
        if len(D) > 0:
            if D[0] == index:
                D.pop(0)
        if len(mis) > 0:
            if mis[0] == index:
                mis.pop(0)
    
    diffs_per_amp = {}
    length = sequences[0].length
    #Deal with lower and upperbound stuff
    if not lb:
        lb = 0
    else:
        lb = max(0,lb)
    if not ub:
        ub = length
    else:
        ub = min(length,ub)
    #Initialize the result dictionary
    diffs_per_amp = dict.fromkeys(set((i,i+amplicon_width) for i in range(lb,ub-amplicon_width+1)))
    for k in diffs_per_amp:
        diffs_per_amp[k] = []
        
    #Initialize list with amplicons to delete due to misaligns
    to_delete = set()
    
    #Iterate over sequence pairs (s1 > s2)
    for s1 in range(len(sequences)):
        for s2 in range(s1):
            D = [] #list containing the indices at which the current pair of sequences differ
            mis = [] #list containing the indices at which either of the pair of sequences contains a misalign
            '''
            Consideration: maybe mis should only be containing the indices at which BOTH sequences contain a misalign?
            currently the situation is less nuanced as both pairs may only have misalign_threshold misalignment characters in aggregation (so not per sequence)
            '''
            #Fill in the first amplicon window
            for i in range(amplicon_width):
                cur_comp = comparison_matrix[(sequences[s1].sequence[i+lb], sequences[s2].sequence[i+lb])]
                if not cur_comp[0]: #cur_comp[0] is True if chars at index are equal, cur_comp[1] is True if either sequence has misalign char '-'
                    D.append(i+lb)
                if cur_comp[1]:
                    mis.append(i+lb)
            i = lb
            #Iterate over the remaining amplicon windows
            while i + amplicon_width < ub:
                #Check for current amplicon whether sequences differ
                if len(mis) > misalign_threshold:
                    to_delete.add((i, i+amplicon_width)) #Amplicons exceeding misalignment criterion will be deleted
                else:
                    #Important to note here, is that if we have the following case (this should be prevented by performing MSA first):
                    #s1 = ~actga-~
                    #s2 = ~-actga~
                    #then s1 and s2 are considered different. The assumption however, is that this is due to an artifact resulting from MSA and shouldn't occur
                    if len(D) > amplicon_threshold:
                        diffs_per_amp[(i,i+amplicon_width)].append((sequences[s2],sequences[s1]))
                    
                #Consider what happens when moving amplicon window 1 index to the right
                moveWindow(D, mis, i, amplicon_width)
                #Check next character
                cur_comp = comparison_matrix[(sequences[s1].sequence[i+amplicon_width], sequences[s2].sequence[i+amplicon_width])]
                if not cur_comp[0]: #If the next character is different, add it to D
                    D.append(i+amplicon_width)
                if cur_comp[1]: #If the next index contains a misalign character, add it to mis
                    mis.append(i+amplicon_width)
                i += 1
            #Consider what happened at the final amplicon window
            moveWindow(D, mis, i, amplicon_width)
            #Check whether misalign condition is satisfied
            if len(mis) > misalign_threshold:
                to_delete.add((i, i+amplicon_width)) #amplicons exceeding misalignment criterion will be deleted
            else:
                if len(D) > amplicon_threshold:
                    diffs_per_amp[(i,i+amplicon_width)].append((sequences[s2],sequences[s1]))
    #Postprocessing: removing all amplicons that contained too many misaligns or unable to differentiate
    res = []
    for amp in diffs_per_amp:
        if not amp in to_delete and len(diffs_per_amp[amp]) > 0:
            temp = Amplicon(amp[0], amp[1])
            temp.setDifferences(diffs_per_amp[amp])
            res.append(temp)
    return res

################################# This code generates and filters primers ##############################################
def generatePrimers(sequences, amplicon, primer_width, search_width, comparison_matrix, at_threshold=3, gc_threshold=4, constrained=True, verbose=False):
    def calcDegeneracy(sequence):
        res = 1
        for char in sequence:
            if char in ['y','r','s','w','m','k']:
                res = res*2
            elif char in ['b','d','h','v']:
                res = res*3
            elif char == 'n':
                res = res*4
        return res
    
    def disambiguate(sequence):
        '''
        Function that "disambiguates" the current (sub)string by generating a string for every degeneracy

        Parameters
        ----------
        sequence : str
            Sequence to "disambiguate".

        Returns
        -------
        res : list[str]
            list containing the disambiguated sequences from the input sequence.

        '''
        translation = {'a' : ['a'],
                       'c' : ['c'],
                       'g' : ['g'],
                       't' : ['t'],                           
                       'b' : ['c','g','t'],
                       'd' : ['a','g','t'],
                       'h' : ['a','c','t'],
                       'k' : ['g','t'],
                       'm' : ['a','c'],
                       'n' : ['a','c','g','t'],
                       'r' : ['a','g'],
                       's' : ['g','c'],
                       'v' : ['a','c','g'],
                       'w' : ['a','t'],
                       'y' : ['c','t']}
        res = translation[sequence[0]].copy()
        for char in sequence[1:]:
            for subsequence_index in range(len(res)):
                new_subsequences = []
                for new_char in translation[char]:
                    new_subsequences.append(res[subsequence_index] + new_char)
                res[subsequence_index] = new_subsequences
            res = list(itertools.chain(*res))
        return res
    
    def reverseComplement(sequence, rev=True):
        #Define complement of every nucleotide
        translate = {
                'a' : 't',
                't' : 'a',
                'u' : 'a',
                'g' : 'c',
                'c' : 'g',
                'y' : 'r',
                'r' : 'y',
                's' : 's',
                'w' : 'w',
                'k' : 'm',
                'm' : 'k',
                'b' : 'v',
                'd' : 'h',
                'h' : 'd',
                'v' : 'b',
                'n' : 'n',
                '-' : '-'
            }
        res = ''
        for i in range(len(sequence)):
            res += translate[sequence[i]] #complement every nucleotide
        if rev:
            return res[::-1]
        else:
            return res
    
    sequence_index = 0
    sequence_amplicon_primers = {}
    forward_primers = {}
    forward_coverage = set()
    reverse_primers = {}
    reverse_coverage = set()
    
    #Updated version:
    for sequence in sequences:
        for offset in range(search_width-primer_width+1):
            #Forward primer:
            fwd_end_index = sequence.aligned_to_trim[amplicon.start-offset-1]
            cur_fwd_seq = sequence.sequence_raw[fwd_end_index-primer_width:fwd_end_index] #get forward primer sequence
            if (calcDegeneracy(cur_fwd_seq) <= 4**6): #Skip if primer is "too degenerate"
                cur_fwd_disambiguated = disambiguate(cur_fwd_seq)
                cur_fwd = []
                for tmp in cur_fwd_disambiguated: #iterate over potential new forward primer candidates
                    if tmp not in forward_primers:
                        forward_primers[tmp] = Primer(tmp, amplicon, 'forward')
            #Reverse primer:
            rev_start_index = sequence.aligned_to_trim[amplicon.end]+offset
            cur_rev_seq = reverseComplement(sequence.sequence_raw[rev_start_index:rev_start_index+primer_width])
            if (calcDegeneracy(cur_rev_seq) <= 4**6): #Skip if primer is "too degenerate"
                cur_rev_disambiguated = disambiguate(cur_rev_seq)
                cur_rev = []
                for tmp in cur_rev_disambiguated:
                    if tmp not in reverse_primers:
                        reverse_primers[tmp] = Primer(tmp, amplicon, 'reverse')
        sequence_index += 1
    if constrained:
        #Filter forward primers
        for primer in forward_primers:
            check = forward_primers[primer].determineProperties(comparison_matrix, at_threshold, gc_threshold)
            if check:
                forward_primers[primer].checkBindingEvents(sequences)
                amplicon.addPrimers(forward_primers[primer], 'forward')
                for sequence in forward_primers[primer].binding_sequences:
                    forward_coverage.add(sequence)
        for primer in reverse_primers:
            check = reverse_primers[primer].determineProperties(comparison_matrix, at_threshold, gc_threshold)
            if check:
                reverse_primers[primer].checkBindingEvents(sequences)
                amplicon.addPrimers(reverse_primers[primer], 'reverse')
                for sequence in reverse_primers[primer].binding_sequences:
                    reverse_coverage.add(sequence)
    else:
        return forward_primers, reverse_primers
    return forward_coverage == set(sequences) and reverse_coverage == set(sequences)

def generatePrimers_exp(sequences, amplicon, primer_width, search_width, comparison_matrix, at_threshold=3, gc_threshold=4, constrained=True, verbose=False):
    def calcDegeneracy(sequence):
        res = 1
        for char in sequence:
            if char in ['y','r','s','w','m','k']:
                res = res*2
            elif char in ['b','d','h','v']:
                res = res*3
            elif char == 'n':
                res = res*4
        return res
    
    def disambiguate(sequence):
        '''
        Function that "disambiguates" the current (sub)string by generating a string for every degeneracy

        Parameters
        ----------
        sequence : str
            Sequence to "disambiguate".

        Returns
        -------
        res : list[str]
            list containing the disambiguated sequences from the input sequence.

        '''
        translation = {'a' : ['a'],
                       'c' : ['c'],
                       'g' : ['g'],
                       't' : ['t'],                           
                       'b' : ['c','g','t'],
                       'd' : ['a','g','t'],
                       'h' : ['a','c','t'],
                       'k' : ['g','t'],
                       'm' : ['a','c'],
                       'n' : ['a','c','g','t'],
                       'r' : ['a','g'],
                       's' : ['g','c'],
                       'v' : ['a','c','g'],
                       'w' : ['a','t'],
                       'y' : ['c','t']}
        res = translation[sequence[0]].copy()
        for char in sequence[1:]:
            for subsequence_index in range(len(res)):
                new_subsequences = []
                for new_char in translation[char]:
                    new_subsequences.append(res[subsequence_index] + new_char)
                res[subsequence_index] = new_subsequences
            res = list(itertools.chain(*res))
        return res
    
    def reverseComplement(sequence, rev=True):
        #Define complement of every nucleotide
        translate = {
                'a' : 't',
                't' : 'a',
                'u' : 'a',
                'g' : 'c',
                'c' : 'g',
                'y' : 'r',
                'r' : 'y',
                's' : 's',
                'w' : 'w',
                'k' : 'm',
                'm' : 'k',
                'b' : 'v',
                'd' : 'h',
                'h' : 'd',
                'v' : 'b',
                'n' : 'n',
                '-' : '-'
            }
        res = ''
        for i in range(len(sequence)):
            res += translate[sequence[i]] #complement every nucleotide
        if rev:
            return res[::-1]
        else:
            return res
    
    sequence_index = 0
    sequence_amplicon_primers = {}
    forward_primers = {}
    forward_coverage = set()
    reverse_primers = {}
    reverse_coverage = set()
    
    #Updated version:
    for sequence in sequences:
        for offset in range(search_width-primer_width):
            #Forward primer:
            fwd_end_index = sequence.aligned_to_trim[amplicon.start-offset]
            cur_fwd_seq = sequence.sequence_raw[fwd_end_index-primer_width:fwd_end_index] #get forward primer sequence
            if (calcDegeneracy(cur_fwd_seq) <= 4**6): #Skip if primer is "too degenerate"
                cur_fwd_disambiguated = disambiguate(cur_fwd_seq)
                cur_fwd = []
                for tmp in cur_fwd_disambiguated: #iterate over potential new forward primer candidates
                    if tmp not in forward_primers and tmp in sequence.kmers['forward']:
                        forward_primers[tmp] = Primer(tmp, amplicon, 'forward')
            #Reverse primer:
            rev_start_index = sequence.aligned_to_trim[amplicon.end]+offset
            cur_rev_seq = reverseComplement(sequence.sequence_raw[rev_start_index:rev_start_index+primer_width])
            if (calcDegeneracy(cur_rev_seq) <= 4**6): #Skip if primer is "too degenerate"
                cur_rev_disambiguated = disambiguate(cur_rev_seq)
                cur_rev = []
                for tmp in cur_rev_disambiguated:
                    if tmp not in reverse_primers and tmp in sequence.kmers['reverse']:
                        reverse_primers[tmp] = Primer(tmp, amplicon, 'reverse')
        sequence_index += 1
    if constrained:
        #Filter forward primers
        for primer in forward_primers:
            forward_primers[primer].checkBindingEvents(sequences)
            amplicon.addPrimers(forward_primers[primer], 'forward')
            for sequence in forward_primers[primer].binding_sequences:
                forward_coverage.add(sequence)
        for primer in reverse_primers:
            reverse_primers[primer].checkBindingEvents(sequences)
            amplicon.addPrimers(reverse_primers[primer], 'reverse')
            for sequence in reverse_primers[primer].binding_sequences:
                reverse_coverage.add(sequence)
    return forward_coverage == set(sequences) and reverse_coverage == set(sequences)

################################# This code contains the greedy helper functions and greedy algorithm #####################
def calcWeight(sequence1, sequence2):
    lin1_split = sequence1.lineage.split('.')
    lin2_split = sequence2.lineage.split('.')
    for i in range(len(lin1_split)):
        if i >= len(lin2_split): #lineage 1 is the same up to i, but longer
            return len(lin1_split) - len(lin2_split)
        elif lin1_split[i] != lin2_split[i]: #lineage 1 is the same up to (not including) i, and different afterwards
            return len(lin1_split[i:]) + len(lin2_split[i:])
    return len(lin2_split) - len(lin1_split) #lineage 2 is the same up to the end of lineage 1, but longer

def calcAmpliconWeight(amplicon, weighted=False):
    res = 0
    if weighted:
        for pair in amplicon.differences:
            res += calcWeight(pair[0], pair[1])
    else:
        res = len(amplicon.differences_proper)
    amplicon.weight = res
    return res

def findMax(amplicons, weighted=False):
    res = (-10**8, None, None)
    for amplicon in amplicons:
        cur_weight = calcAmpliconWeight(amplicon, weighted=weighted)
        if cur_weight > res[0] and len(amplicon.differences) > 0:
            res = (cur_weight, amplicon, amplicon.differences.copy())
    return res

#Deprecated version of greedyWithPrimers
'''
def greedyWithPrimers(sequences, amplicons, primer_width, search_width, comparison_matrix, max_iterations=100, weighted=False, logging=None):
    to_cover = set()
    for amplicon in amplicons:
        amplicon.differences_proper = set()
        for pair in amplicon.differences:
            if pair[0].lineage != pair[1].lineage:
                amplicon.differences_proper.add(pair)
                to_cover.add(pair)
    total_to_cover = len(to_cover)
    if logging:
        log_results = ['Total to cover: ' + str(total_to_cover) + ' out of ' + str( len(sequences) * (len(sequences)-1) / 2 )]
    
    results_amplicons = []
    results_primers = []
    remaining = copy.deepcopy(amplicons)
    while len(to_cover) > 0 and len(results_amplicons) < max_iterations:
        cur_max = findMax(remaining, weighted=weighted)
        to_remove = []
        print('Checking amplicon: ' + str(cur_max[1].id))
        if logging:
            log_results.append('Checking amplicon: ' + str(cur_max[1].id))
        if generatePrimers(sequences, cur_max[1], primer_width, search_width, comparison_matrix):
            results_amplicons.append(cur_max[1])
            cur_primers = checkFeasibility(results_amplicons, sequences, comparison_matrix, optimize=0, temperature_range=5)
            if cur_primers:
                to_cover = to_cover - cur_max[1].differences
                remaining.remove(cur_max[1])
                results_primers.append(cur_primers)
                print('Accepted: ' + str(cur_max[1].id))
                if logging:
                    log_results.append('Amplicon ' + str(cur_max[1].id) + ' succesfully added, coverage=' + str(len(cur_max[1].differences)))
                for amplicon in remaining:
                    amplicon.differences = amplicon.differences - cur_max[1].differences
                    amplicon.differences_proper = amplicon.differences_proper - cur_max[1].differences
                    if len(amplicon.differences) == 0 or abs(amplicon.start - cur_max[1].start) < (amplicon.end - amplicon.start) + search_width:
                        to_remove.append(amplicon)
            else:
                results_amplicons.pop(-1)
                print('Rejected: ' + str(cur_max[1].id))
                if logging:
                    log_results.append('Amplicon ' + str(cur_max[1].id) + ' rejected as no feasible primers can be found')
                to_remove.append(cur_max[1])
            for amplicon in to_remove:
                remaining.remove(amplicon)
        else:
            print('Unable to find primers for amplicon: ' + str(cur_max[1].id))
            if logging:
                log_results.append('Amplicon ' + str(cur_max[1].id) + ' rejected as primers do not cover every sequence')
            remaining.remove(cur_max[1])
    if logging:
        return log_results
    else:
        return results_amplicons, results_primers
'''
    
def greedyWithPrimers_exp(sequences, amplicons, primer_width, search_width, comparison_matrix, max_iterations=100, weighted=False, logging=None, multiplex=True):
    to_cover = set()
    for amplicon in amplicons:
        amplicon.differences_proper = set()
        for pair in amplicon.differences:
            if pair[0].lineage != pair[1].lineage:
                amplicon.differences_proper.add(pair)
                to_cover.add(pair)
    total_to_cover = len(to_cover)
    if logging:
        log_results = ['Total to cover: ' + str(total_to_cover) + ' out of ' + str( len(sequences) * (len(sequences)-1) / 2 )]
    results_amplicons = []
    amplicons.sort(key = lambda x: len(x.differences_proper), reverse=True) #Sort in descending order in terms of correct differences
    removed_amplicons = [] #Stores removed amplicons so that they can be added back later
    while len(to_cover) > 0 and len(results_amplicons) < max_iterations and len(amplicons) > 0:
        cur_best_amplicon = amplicons.pop(0) #Remove current best from list
        removed_amplicons.append(cur_best_amplicon) #Move to list of removed amplicons
        to_remove = []
        #print('Checking amplicon: ' + str(cur_best_amplicon.id))
        if logging:
            log_results.append('Checking amplicon: ' + str(cur_best_amplicon.id))
            print('Checking amplicon: ' + str(cur_best_amplicon.id))
        if generatePrimers(sequences, cur_best_amplicon, primer_width, search_width, comparison_matrix):
            results_amplicons.append(cur_best_amplicon)
            if multiplex:
                cur_primer_check, _ = checkFeasibility(results_amplicons, sequences, comparison_matrix, optimize=0, temperature_range=5)
            else:
                cur_primer_check, _ = checkFeasibility([cur_best_amplicon], sequences, comparison_matrix, optimize=0, temperature_range=5)
            if cur_primer_check:
                to_cover = to_cover - cur_best_amplicon.differences
                #print('Accepted: ' + str(cur_best_amplicon.id))
                if logging:
                    log_results.append('Amplicon ' + str(cur_best_amplicon.id) + ' succesfully added, coverage=' + str(len(cur_best_amplicon.differences_proper)))
                for amplicon in amplicons:
                    amplicon.differences_proper = amplicon.differences_proper - cur_best_amplicon.differences
                    if len(amplicon.differences_proper) == 0 or abs(amplicon.start - cur_best_amplicon.start) < (amplicon.end - amplicon.start) + search_width:
                        to_remove.append(amplicon)
                for amplicon in to_remove:
                    removed_amplicons.append(amplicon)
                    amplicons.remove(amplicon)
                amplicons.sort(key = lambda x: len(x.differences_proper), reverse=True)
            else:
                results_amplicons.pop(-1)
                #print('Rejected: ' + str(cur_best_amplicon.id))
                if logging:
                    log_results.append('Amplicon ' + str(cur_best_amplicon.id) + ' rejected as no feasible primers can be found')
        else:
            #print('Unable to find primers for amplicon: ' + str(cur_best_amplicon.id))
            if logging:
                log_results.append('Amplicon ' + str(cur_best_amplicon.id) + ' rejected as primers do not cover every sequence')
    amplicons = amplicons + removed_amplicons
    if logging:
        return log_results, amplicons, results_amplicons
    else:
        return results_amplicons, results_primers

def checkFeasibility(amplicons, sequences, comparison_matrix, optimize=0, temperature_range=5):
    env = gp.Env(empty=True)
    env.setParam('OutputFlag',0) #Turn off logging in interactive window
    env.start()
    
    model = gp.Model(env=env)
    model.ModelSense=GRB.MINIMIZE
    forward_primers = {} #Stores values of forward primers, primer.sequence -> primer variable in ILP
    forward_temperatures = {}
    reverse_primers = {} #Stores values of reverse primers, primer.sequence -> primer variable in ILP
    reverse_temperatures = {}
    
    #Initialize Gurobi variables
    for amplicon in amplicons:
        for primer in amplicon.primers['forward']:
            if primer.sequence not in forward_primers:
                forward_primers[primer.sequence] = model.addVar(vtype=GRB.BINARY, obj=optimize)
                forward_temperatures[primer.sequence] = primer.melting_temperature #Store melting temperature
        for primer in amplicon.primers['reverse']:
            if primer.sequence not in reverse_primers:
                reverse_primers[primer.sequence] = model.addVar(vtype=GRB.BINARY, obj=optimize)
                reverse_temperatures[primer.sequence] = primer.melting_temperature #Store melting temperature
    
    #Every sequence should be covered by both a forward and a reverse primer for every amplicon
    #Note that checkFeasibility is only run when we have an amplicon that can be amplified in every sequence
    for amplicon in amplicons:
        for sequence in amplicon.primers_per_sequence:
            model.addConstr(sum(forward_primers[primer.sequence] for primer in amplicon.primers_per_sequence[sequence]['forward']) >= 1)
            model.addConstr(sum(reverse_primers[primer.sequence] for primer in amplicon.primers_per_sequence[sequence]['reverse']) >= 1)
            
    #Managing conflicting primers
    full_primer_list = set()
    for amplicon in amplicons:
        for primer in amplicon.primers['forward']:
            full_primer_list.add(primer)
        for primer in amplicon.primers['reverse']:
            full_primer_list.add(primer)
    #Iterate through every primer pair and check if they conflict
    for primer_pair in itertools.combinations(full_primer_list, 2):
        if primer_pair[0].checkCompatibility(primer_pair[1], comparison_matrix)[0] >= 11:
            if primer_pair[0].orientation == 'forward' and primer_pair[1].orientation == 'forward':
                model.addConstr(forward_primers[primer_pair[0].sequence] + forward_primers[primer_pair[1].sequence] <= 1)
            elif primer_pair[0].orientation == 'forward' and primer_pair[1].orientation == 'reverse':
                model.addConstr(forward_primers[primer_pair[0].sequence] + reverse_primers[primer_pair[1].sequence] <= 1)
            elif primer_pair[0].orientation == 'reverse' and primer_pair[1].orientation == 'forward':
                model.addConstr(reverse_primers[primer_pair[0].sequence] + forward_primers[primer_pair[1].sequence] <= 1)
            else:
                model.addConstr(reverse_primers[primer_pair[0].sequence] + reverse_primers[primer_pair[1].sequence] <= 1)
                
    #Enforce that the max melting temperature and min melting temperature are within specified range
    min_temp = model.addVar(vtype=GRB.CONTINUOUS)
    max_temp = model.addVar(vtype=GRB.CONTINUOUS)
    for primer in forward_primers:
        model.addConstr(min_temp <= forward_temperatures[primer] * (3 - 2*forward_primers[primer]))
        model.addConstr(max_temp >= forward_temperatures[primer] * forward_primers[primer])
    for primer in reverse_primers:
        model.addConstr(min_temp <= reverse_temperatures[primer] * (3 - 2*reverse_primers[primer]))
        model.addConstr(max_temp >= reverse_temperatures[primer] * reverse_primers[primer])
    model.addConstr(max_temp - min_temp <= temperature_range)
    
    model.optimize()
    if model.Status == 2 and optimize > 0:
        final_primers = {'forward' : [], 'reverse' : []}
        for primer in forward_primers:
            if forward_primers[primer].x >= 0.01:
                final_primers['forward'].append(primer)
        for primer in reverse_primers:
            if reverse_primers[primer].x >= 0.01:
                final_primers['reverse'].append(primer)
        return (model.Status == 2), final_primers
    
    return model.Status == 2, None

def checkFeasibility_naive_percentage(amplicons, sequences, comparison_matrix, optimize=0, temperature_range=5, coverage=0):
    env = gp.Env(empty=True)
    env.setParam('OutputFlag',0) #Turn off logging in interactive window
    env.start()
    
    model = gp.Model(env=env)
    model.ModelSense=GRB.MINIMIZE
    forward_primers = {} #Stores values of forward primers, primer.sequence -> primer variable in ILP
    forward_temperatures = {}
    reverse_primers = {} #Stores values of reverse primers, primer.sequence -> primer variable in ILP
    reverse_temperatures = {}
    binary_check = {} #Stores values representing whether a sequence amplicon pair is covered
    
    #Initialize Gurobi variables
    for amplicon in amplicons:
        for primer in amplicon.primers['forward']:
            if primer.sequence not in forward_primers:
                forward_primers[primer.sequence] = model.addVar(vtype=GRB.BINARY, obj=optimize)
                forward_temperatures[primer.sequence] = primer.melting_temperature #Store melting temperature
        for primer in amplicon.primers['reverse']:
            if primer.sequence not in reverse_primers:
                reverse_primers[primer.sequence] = model.addVar(vtype=GRB.BINARY, obj=optimize)
                reverse_temperatures[primer.sequence] = primer.melting_temperature #Store melting temperature
        for sequence in sequences:
            binary_check[(amplicon.id, sequence.id)] = model.addVar(vtype=GRB.BINARY, obj=0)
    
    #Every sequence should be covered by both a forward and a reverse primer for every amplicon
    #Note that checkFeasibility is only run when we have an amplicon that has potential primers for every sequence
    for amplicon in amplicons:
        for sequence in amplicon.primers_per_sequence:
            model.addConstr(binary_check[(amplicon.id, sequence.id)] <= sum(forward_primers[primer.sequence] for primer in amplicon.primers_per_sequence[sequence]['forward']))
            model.addConstr(binary_check[(amplicon.id, sequence.id)] <= sum(reverse_primers[primer.sequence] for primer in amplicon.primers_per_sequence[sequence]['reverse']))
        model.addConstr(sum(binary_check[(amplicon.id, sequence.id)] for sequence in amplicon.primers_per_sequence) >= coverage * len(sequences))
        
    #Managing conflicting primers
    full_primer_list = set()
    for amplicon in amplicons:
        for primer in amplicon.primers['forward']:
            full_primer_list.add(primer)
        for primer in amplicon.primers['reverse']:
            full_primer_list.add(primer)
    #Iterate through every primer pair and check if they conflict
    for primer_pair in itertools.combinations(full_primer_list, 2):
        if primer_pair[0].checkCompatibility(primer_pair[1], comparison_matrix)[0] >= 9:
            if primer_pair[0].orientation == 'forward' and primer_pair[1].orientation == 'forward':
                model.addConstr(forward_primers[primer_pair[0].sequence] + forward_primers[primer_pair[1].sequence] <= 1)
            elif primer_pair[0].orientation == 'forward' and primer_pair[1].orientation == 'reverse':
                model.addConstr(forward_primers[primer_pair[0].sequence] + reverse_primers[primer_pair[1].sequence] <= 1)
            elif primer_pair[0].orientation == 'reverse' and primer_pair[1].orientation == 'forward':
                model.addConstr(reverse_primers[primer_pair[0].sequence] + forward_primers[primer_pair[1].sequence] <= 1)
            else:
                model.addConstr(reverse_primers[primer_pair[0].sequence] + reverse_primers[primer_pair[1].sequence] <= 1)
                
    #Enforce that the max melting temperature and min melting temperature are within specified range
    min_temp = model.addVar(vtype=GRB.CONTINUOUS)
    max_temp = model.addVar(vtype=GRB.CONTINUOUS)
    for primer in forward_primers:
        model.addConstr(min_temp <= forward_temperatures[primer] * (3 - 2*forward_primers[primer]))
        model.addConstr(max_temp >= forward_temperatures[primer] * forward_primers[primer])
    for primer in reverse_primers:
        model.addConstr(min_temp <= reverse_temperatures[primer] * (3 - 2*reverse_primers[primer]))
        model.addConstr(max_temp >= reverse_temperatures[primer] * reverse_primers[primer])
    model.addConstr(max_temp - min_temp <= temperature_range)
    
    model.optimize()
    if model.Status == 2 and optimize > 0:
        final_primers = {'forward' : [], 'reverse' : []}
        for primer in forward_primers:
            if forward_primers[primer].x >= 0.01:
                final_primers['forward'].append(primer)
        for primer in reverse_primers:
            if reverse_primers[primer].x >= 0.01:
                final_primers['reverse'].append(primer)
        return (model.Status == 2), final_primers
    
    return model.Status == 2, None
################################# This code deals with visualization and stuff ############################################
"""
    
    
    
    
    
    
    
    
    