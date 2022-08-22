import RNA
from Bio.SeqUtils import MeltingTemp as mt
import itertools
import os
import multiprocess as mp
import numpy as np
from Scripts import *
from math import ceil
import argparse

"""
To-do:
    - add to_text methods to custom classes
    - add read methods to custom classes <- not sure if this is necessary
    - update generate_primers
    - update generate_amplicons
    - update greedy
"""

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


#Updated classes
class Sequence:
    def __init__(self, sequence, identifier, lineage=None):
        self.sequence = sequence
        self.length = len(sequence)
        self.sequence_raw = sequence.replace('-','')
        self.length_raw = len(self.sequence_raw)
        self.id = identifier
        self.lineage = lineage
        self.aligned_to_trim = np.zeros((1))
        
    def __eq__(self, other):
        try:
            if type(other) == str:
                return self.sequence == other or self.sequence_raw == other or self.id == other
            else:
                return self.id == other.id
        except:
            return False
        
    def __hash__(self):
        return hash(self.id)
    
    def __repr__(self):
        return 'Sequence'
    
    def align_to_trim(self):
        '''
        Function that determines which character index in the original sequence refers to which index in the raw sequence

        Returns
        -------
        np.array[int]
            np.array with indices representing the index of a character in the original sequence, and the element is equal to the corresponding element in the raw sequence.

        '''
        self.aligned_to_trim = np.zeros((self.length), dtype=int)
        map_index = -1 #note that there is always a problem with how we define the "first" index or the "last" index
        for char_index in range(self.length):
            if self.sequence[char_index] != '-':
                map_index += 1
            self.aligned_to_trim[char_index] = map_index
        return self.aligned_to_trim
    
    def find_bounds(self, min_non_align):
        '''
        Function that determines the number of "nucleotides" that have to preceed (exceed) such that there are at exactly $min_non_align nucleotides before (after) it. This requires
        the existence of the aligned_to_trim list in the Sequence object in order to work.
        
        As an example consider the sequence --ac-t-gtacctg-a-gc- ([-1, -1, 0, 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9, 10, 10, 11, 12, 12]) and min_non_align=3.
        Since we want at least 3 non-align characters, we should find the first index in the aligned_to_trim list equal to 3 as before this index there are exactly 3 non-align characters.
        Similarly, there is a total of 13 non-align chars, and we thus need the index of the 13-3=10th non-align such that everything afterwards contains exactly 3 non-align characters.
        It should be noted that there may be some ambiguity in this because as in the example there are two elements equal to 10. The goal is to find lb and ub such that:
            self.sequence[:lb] and self.sequence[ub:] both contain exactly &min_non_align non-align characters and,
            moving lb one to the left, or ub one to the right no longer satisfies the previous property

        Parameters
        ----------
        min_non_align : int
            Number of non-alignment characters before and after the bounds to be found.

        Returns
        -------
        lb : int
            index such that there are exactly $min_non_align non-alignment characters before it.
        ub : int
            index such that there are exactly $min_non_align non-alignment characters before it.

        '''
        #lb = self.aligned_to_trim.index(min_non_align) <- old
        lb = np.where(self.aligned_to_trim == min_non_align)[0][0]
        #ub = self.aligned_to_trim.index(self.length_raw - min_non_align) <- old
        ub = np.where(self.aligned_to_trim == (self.length_raw - min_non_align))[0][0]
        return lb, ub

class Amplicon:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.id = (start,end)
        #Note here that primers are stored as their index in the PrimerIndex, not as Primer object, and every orientation holds the primer indices per sequence
        self.primers = {'forward' : {}, 'reverse': {}} 
        self.differences = set()
        self.differences_proper = set()
        
    def __eq__(self, other):
        try:
            if type(other) == tuple:
                return self.id == other
            else:
                return self.id == other.id
        except:
            return False
        
    def __lt__(self, other):
        #Note here that we can either check for the amplicon index, or the differentiation power. Here we use the latter
        return len(self.differences_proper) < len(other.differences_proper)
    
    def __repr__(self):
        return 'Amplicon'
    
    def add_primers(self, primers):
        '''
        Function that adds primers (as strings) to current Amplicon object

        Parameters
        ----------
        primers : dict[ str ][ str ]
            Dictionary with orientation as primary key, sequence ids as secondary key and the primer indices as values

        Returns
        -------
        None.

        '''
        for orientation in primers:
            for sequence in primers[orientation]:
                if sequence in self.primers[orientation]:
                    self.primers[orientation][sequence] = self.primers[orientation][sequence].union(primers[orientation][sequence])
                else:
                    self.primers[orientation][sequence] = primers[orientation][sequence]
        
    def set_differences(self, differences):
        '''
        Function that sets the differences of this Amplicon object to the differences given. Additionally, it tries to do the same but only for "proper" differences.

        Parameters
        ----------
        differences : set[Sequence], list[Sequence]
            Set or list of pairs of Sequences objects indicating which sequence pairs can be differentiated based on this amplicon.

        Returns
        -------
        None.

        '''
        self.differences = set(differences)
        try:
            for difference in differences:
                if difference[0].lineage != difference[1].lineage:
                    self.differences_proper.add((difference[0].id, difference[1].id))
        except:
            pass
        
    def check_differences(self, sequences):
        for difference in self.differences:
            try:
                if sequences[sequences.index(difference[0])].lineage != sequences[sequences.index(difference[1])].lineage:
                    self.differences_proper.add((difference[0], difference[1]))
            except:
                continue
       
    #Deprecated function
    '''
    def set_weight(self, lineages, ignore, weighted=False):
        def calc_weight(lineage1, lineage2):
            lin1_split = lineage1.split('.')
            lin2_split = lineage2.split('.')
            for i in range(len(lin1_split)):
                if i >= len(lin2_split):
                    return len(lin1_split) - len(lin2_split)
                elif lin1_split[i] != lin2_split[i]:
                    return len(lin1_split[i:]) + len(lin2_split[i:])
            return len(lin2_split) - len(lin1_split)
        
        self.weight = 0
        if weighted:
            for pair in (self.differences - ignore):
                self.weight += calcWeight(lineages[pair[0]], lineages[pair[1]])
        else:
            self.weight = len(self.differences - ignore)
        return self.weight
    '''
                
class Primer:
    def __init__(self, sequence, orientation):
        self.sequence = sequence
        self.indices = {} #stores the starting indices of this primer for the sequences of interest
        self.orientation = orientation
        self.feasible = True #By default, feasibility of a primer is set to true
        self.temperature = None
        
    def __eq__(self, other):
        try:
            if type(other) == str:
                return other == self.sequence
            else:
                return other.sequence == self.sequence
        except:
            return False
        
    def add_sequence(self, sequence, sequence_index):
        '''
        Function that adds the occurrence of this primer to a particular sequence object

        Parameters
        ----------
        sequence : Sequence
            Sequence object that contains this primer.
        sequence_index : int
            Starting index of the primer in sequence (e.g. primer=act, sequence=ggactca -> starting_index=2 as sequence[2:2+3] = act).

        Returns
        -------
        bool
            True if the sequence has been succesfully added, False otherwise.

        '''
        if self.feasible:
            if sequence.id not in self.indices:
                self.indices[sequence.id] = sequence_index
            else:
                if not self.indices[sequence.id] == sequence_index:
                    self.indices = None
                    self.feasible = False
                    return False
        return True
            
    def check_compatibility(self, other, comparison_matrix, comparison_tolerance):
        '''
        Function that checks whether this primer is compatible with the $other primer based on complementarity in all possible alignments.

        Parameters
        ----------
        other : Primer
            Primer object to check compatibility with.
        comparison_matrix : dict[ (char,char) ]
            Dictionary that determines which characters should be considered equal.
        comparison_tolerance : int
            Maximal allowed complementarity between self and $other.

        Returns
        -------
        res : (int, str)
            Tuple where the first element is the maximum complementarity (or the first complementarity that exceeds the tolerance), and
            the second element is a visual representation of this complementarity.

        '''
        res = (0, '')
        to_compare = reverse_complement(other.sequence)
        
        for i in range(len(self.sequence)):
            fwd_check = []
            rev_check = []
            for j in range(i, len(self.sequence)):
                #Check forward comparison
                if comparison_matrix[(to_compare[j], self.sequence[j-i])][0] and not comparison_matrix[(to_compare[j], self.sequence[j-i])][1]:
                    fwd_check.append(1)
                else:
                    fwd_check.append(0)
                
                #Check backward comparison
                if comparison_matrix[(self.sequence[j], to_compare[j-i])][0] and not comparison_matrix[(self.sequence[j], to_compare[j-i])][1]:
                    rev_check.append(1)
                else:
                    rev_check.append(0)
            #Check forward
            comparison = sum(fwd_check)
            if comparison > res[0]:
                if self == other:
                    res = [comparison, i*'*' + self.sequence + '\n' + self.sequence[::-1] + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in fwd_check) + i*'*']
                else:
                    res = [comparison, i*'*' + self.sequence + '\n' + other.sequence[::-1] + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in fwd_check) + i*'*']
            #If number of matches exceeds upperbound then return
            if comparison > comparison_tolerance:
                return res
            #Check reverse
            comparison = sum(rev_check)
            if comparison > res[0]:
                if self == other:
                    res = [comparison, i*'*' + self.sequence[::-1] + '\n' + self.sequence + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in rev_check) + i*'*']
                else:
                    res = [comparison, i*'*' + other.sequence[::-1] + '\n' + self.sequence + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in rev_check) + i*'*']
            #If number of matches exceeds upperbound then return
            if comparison > comparison_tolerance:
                return res
        return res
            
    def check_feasibility(self, comparison_matrix, gc_lb=0.4, gc_ub=0.6, melting_lb=55, melting_ub=75, end_at_threshold=2, end_gc_threshold=3, monorun_threshold=3, duorun_threshold=2, mfe_threshold=-5, self_complementarity_threshold=10, verbose=False):
        '''
        Function that determines whether this primer object is feasible given several chemical and thermal constraints.

        Parameters
        ----------
        comparison_matrix : dict[ (char,char) ]
            Dictionary that determines which characters should be considered equal.
        gc_lb : float, optional
            Minimal GC-content. The default is 0.4.
        gc_ub : float, optional
            Maximal GC-content. The default is 0.6.
        melting_lb : float, optional
            Minimal melting temperature. The default is 55.
        melting_ub : float, optional
            Maximal melting temperature. The default is 75.
        end_at_threshold : int, optional
            Maximal number of allowed A/T characters in final 3 nucleotides (3'-end). The default is 2.
        end_gc_threshold : int, optional
            Maximal number of allowed C/G characters in final 5 nucleotides (5'-end). The default is 3.
        monorun_threshold : int, optional
            Maximal length of a run of a single nucleotide. The default is 3.
        duorun_threshold : int, optional
            Maximal length of a run of a pair of nucleotides. The default is 2.
        mfe_threshold : float, optional
            Minimal allowed MFE representing the "risk" of hairpin formation. The default is -5.
        self_complementarity_threshold : int, optional
            Maximal number of complementary base-pairs in every alignment of this primer with its own reverse. The default is 10.
        verbose : bobl, optional
            True if everything should be saved and outputted, False for fast check. The default is False.

        Returns
        -------
        bool
            True if this primer is feasible by itself, False otherwise.

        '''
        if verbose:
            problems = []
            self.degeneracy = calculate_degeneracy(self.sequence)
            if self.degeneracy > 1:
                problems.append('Primer is too degenerate')
            self.gc = calculate_GC(self.sequence)
            if self.gc < gc_lb:
                problems.append('GC content is too low')
            elif self.gc > gc_ub:
                problems.append('GC content is too high')
            self.end_stats = calculate_end_stats(self.sequence, comparison_matrix)
            if self.end_stats[0] > end_at_threshold:
                problems.append('Too many a/t in final 3 nucleotides in 3\' end')
            if self.end_stats[1] > end_gc_threshold:
                problems.append('Too many g/c in final 5 nucleotides in 5\' end')
            self.mono_run = calculate_longest_monorun(self.sequence, comparison_matrix)
            if self.mono_run > monorun_threshold:
                problems.append('Longest monorun exceeds threshold')
            self.duo_run = calculate_longest_duorun(self.sequence, comparison_matrix)
            if self.duo_run > duorun_threshold:
                problems.append('Longest duorun exceeds threshold')
            try:
                self.temperature = mt.Tm_NN(self.sequence, strict=False)
            except:
                self.temperature = 10000
            if self.temperature < melting_lb:
                problems.append('Melting temperature is too low')
            elif self.temperature > melting_ub:
                problems.append('Melting temperature is too high')
            self.self_complementarity = self.check_compatibility(self, comparison_matrix, self_complementarity_threshold)
            if self.self_complementarity[0] > self_complementarity_threshold:
                problems.append('Self complementarity risk too high')
            #Print all issues with the current primer and then return whether it is feasible
            for p in problems: print(p)
            return len(problems) == 0
        else:
            if calculate_degeneracy(self.sequence) > 1:
                self.feasible = False
                return False
            cur_crit = calculate_GC(self.sequence)
            if cur_crit < gc_lb or cur_crit > gc_ub:
                self.feasible = False
                return False
            cur_crit = calculate_end_stats(self.sequence, comparison_matrix)
            if cur_crit[0] > end_at_threshold or cur_crit[1] > end_gc_threshold:
                self.feasible = False
                return False
            if calculate_longest_monorun(self.sequence, comparison_matrix) > monorun_threshold:
                self.feasible = False
                return False
            if calculate_longest_duorun(self.sequence, comparison_matrix) > duorun_threshold:
                self.feasible = False
                return False
            self.temperature = mt.Tm_NN(self.sequence, strict=False)
            if self.temperature < melting_lb or self.temperature > melting_ub:
                self.feasible = False
                return False
            if self.check_compatibility(self, comparison_matrix, self_complementarity_threshold)[0] > self_complementarity_threshold:
                self.feasible = False
                return False
            return True
            
class PrimerIndex():
    def __init__(self):
        self.primer2index = {'forward' : {}, 'reverse' : {}} #contains primer sequences as keys and corresponding primer index as value
        self.index2primer = {'forward' : np.empty((0)), 'reverse' : np.empty((0))} #contains Primer objects in an array
        self.conflict_matrix = None
        
        self.thresholds = {
            'gc_lb'                             : 0.4,
            'gc_ub'                             : 0.6,
            'melting_lb'                        : 55,
            'melting_ub'                        : 75,
            'end_at_threshold'                  : 2,
            'end_gc_threshold'                  : 3,
            'monorun_threshold'                 : 3,
            'duorun_threshold'                  : 3,
            'mfe_threshold'                     : -5,
            'self_complementarity_threshold'    : 10
            }
        self.comparison_matrix = generate_opportunistic_matrix()
        
    #This does not work currently
    def __eq__(self, other):
        try:
            for orientation in self.set:
                for p in self.set[orientation]:
                    if not self.set[orientation][p].indices == other.set[orientation][p].indices or self.set[orientation][p].feasible == other.set[orientation][p].feasible:
                        return False
            return True
        
            for orientation in other.set:
                for p in other.set[orientation]:
                    if not self.set[orientation][p].indices == other.set[orientation][p].indices or self.set[orientation][p].feasible == other.set[orientation][p].feasible:
                        return False
            return True
        except:
            return False
        
    #Remove new if confirmed to work
    def add_primer(self, sequence, orientation):
        '''
        Function that adds a new primer to the primer index

        Parameters
        ----------
        sequence : str
            String representation of the primer to add.
        orientation : str
            'forward' or 'reverse' indicating the primer orientation.

        Returns
        -------
        None.

        '''
        if sequence not in self.primer2index[orientation]:
            self.primer2index[orientation][sequence] = len(self.index2primer[orientation])
            self.index2primer[orientation] = np.append(self.index2primer[orientation], Primer(sequence, orientation))
            #self.index2primer[orientation].append(Primer(sequence, orientation))
            self.index2primer[orientation][-1].check_feasibility(self.comparison_matrix,
                                                              gc_lb=self.thresholds['gc_lb'],
                                                              gc_ub=self.thresholds['gc_ub'],
                                                              melting_lb=self.thresholds['melting_lb'],
                                                              melting_ub=self.thresholds['melting_ub'],
                                                              end_at_threshold=self.thresholds['end_at_threshold'],
                                                              end_gc_threshold=self.thresholds['end_gc_threshold'],
                                                              monorun_threshold=self.thresholds['monorun_threshold'],
                                                              duorun_threshold=self.thresholds['duorun_threshold'],
                                                              mfe_threshold=self.thresholds['mfe_threshold'],
                                                              self_complementarity_threshold=self.thresholds['self_complementarity_threshold'])

    def add_sequence(self, sequence, sequence_index, primer_sequence, orientation):
        '''
        Function that adds a sequence to a primer in the index, or create a new primer and add the sequence to it

        Parameters
        ----------
        sequence : Sequence
            Sequence object that should be related to a primer.
        sequence_index : str
            Starting index of the primer that has to be linked to this sequence.
        primer_sequence : str
            String representation of the primer to add.
        orientation : str
            'forward' or 'reverse' indicating the primer orientation.

        Returns
        -------
        bool
            True if the sequence has been succesfully added, False otherwise

        '''
        if primer_sequence not in self.primer2index[orientation]:
            self.add_primer(primer_sequence, orientation)
            self.index2primer[orientation][-1].add_sequence(sequence, sequence_index)
        else:
            self.index2primer[orientation][self.primer2index[orientation][primer_sequence]].add_sequence(sequence, sequence_index)
        return self.index2primer[orientation][-1].feasible

    def remove_redundant(self):
        '''
        Function that removes all of the primers in this PrimerIndex that are infeasible. Note that the function itself
        does not check feasibility, instead see Primer.check_feasibility!

        Returns
        -------
        None.

        '''
        common_kmers = set(self.primer2index['forward'].keys()).intersection(set(self.primer2index['reverse'].keys()))
        for orientation in self.primer2index:
            kmers = list(self.primer2index[orientation].keys())
            print( 'Initially contains %d %s primers' % (len(kmers), orientation) )
            index = 0
            to_remove = []
            while index < len(kmers):
                if not self.index2primer[orientation][index].feasible or kmers[index] in common_kmers:
                    to_remove.append(index)
                    self.primer2index[orientation].pop(kmers[index])
                index += 1
            self.index2primer[orientation] = np.delete(self.index2primer[orientation], to_remove)
            for index in range(len(self.index2primer[orientation])):
                self.primer2index[orientation][self.index2primer[orientation][index].sequence] = index
            print( 'Finally contains %d %s primers' % (len(self.primer2index[orientation]), orientation) )
            print( 'Removed %d primers occurring both as forward and reverse' % (len(common_kmers)) )
   
    def set_thresholds(self, thresholds):
        '''
        Function that sets the primer property thresholds to the given values. Note that this does not check for
        existing primers in the index whether they satisfy the new thresholds.

        Parameters
        ----------
        thresholds : dict[ String ]
            Dictionary containing the properties as keys, and the values they should be set to as values.

        Returns
        -------
        None.

        '''
        for prop in thresholds:
            try:
                self.thresholds[prop] = thresholds[prop]
            except:
                continue

    def merge_indices(self, other_index):
        for orientation in other_index.primer2index:
            primers_to_add = []
            indices_to_add = []
            for primer in other_index.primer2index[orientation]:
                index_in_other = other_index.primer2index[orientation][primer]
                #Check if primer is also in this index
                if primer in self.primer2index[orientation]:
                    index_in_this = self.primer2index[orientation][primer]
                    #Check if primer is feasible both in this index and other index
                    if self.index2primer[orientation][index_in_this].feasible and other_index.index2primer[orientation][index_in_other].feasible:
                        for sequence in other_index.index2primer[orientation][index_in_other].indices:
                            #Check if the sequence is also linked to the primer in this index
                            if sequence in self.index2primer[orientation][self.primer2index[orientation][primer]].indices:
                                #Check if primer occurs at the same indices, otherwise it should be rejected
                                if self.index2primer[orientation][index_in_this].indices[sequence] == other_index.index2primer[orientation][index_in_other].indices[sequence]:
                                    continue
                                else:
                                    self.index2primer[orientation][index_in_this].feasible = False
                            #If sequence is not linked to the primer in this, add it
                            else:
                                self.index2primer[orientation][index_in_this].indices[sequence] = copy.deepcopy(other_index.index2primer[orientation][index_in_other].indices[sequence])
                    #If primer is infeasible in either index then set to infeasible
                    else:
                        self.index2primer[orientation][index_in_this].feasible = False
                #If primer is not in this index add it and copy information from other index
                else:
                    primers_to_add.append(primer)
                    indices_to_add.append(index_in_other)
            k = 0
            for primer in primers_to_add:
                self.primer2index[orientation][primer] = len(self.index2primer[orientation]) + k
                k += 1
            self.index2primer[orientation] = np.append(self.index2primer[orientation], other_index.index2primer[orientation][indices_to_add])

    def check_amplicon(self, sequences, amplicon, primer_width, search_width):
        '''
        Function that generates the primers (per sequence) of length $primer_width in a search window of length $search_width
        that can be used to amplify this amplicon for all the sequences in $sequences.

        Parameters
        ----------
        sequences : list[ Sequence ]
            List of sequences to generate primers for.
        amplicon : Amplicon
            Amplicon to find primers around.
        primer_width : int
            Width of primers in number of nucleotides.
        search_width : int
            Search window around the amplicon in which we want to find primers.

        Returns
        -------
        None
        
        '''
        sequence_ids = [sequence.id for sequence in sequences]
        amplicon.primers = {'forward' : {s : set() for s in sequence_ids}, 'reverse' : {s : set() for s in sequence_ids}}
        amplicon.full_primerset = {'forward' : set(), 'reverse' : set()}
        for sequence in sequences:
            #Check if the start of the amplicon is a misalign, in which case correct for it
            if sequence.aligned_to_trim[amplicon.start] == sequence.aligned_to_trim[amplicon.start-1]:
                forward_end_index = sequence.aligned_to_trim[amplicon.start] + 1
            else:
                forward_end_index = sequence.aligned_to_trim[amplicon.start]
            #Check if the character after the amplicon is a misalign, in which case correct for it
            if sequence.aligned_to_trim[amplicon.end] == sequence.aligned_to_trim[amplicon.end-1]:
                reverse_start_index = sequence.aligned_to_trim[amplicon.end] + 1
            else:
                reverse_start_index = sequence.aligned_to_trim[amplicon.end]
                
            #Check which primers (both forward and reverse) are found for the corresponding sequences
            for offset in range(search_width - primer_width + 1): #this iterates over the possible primers within the search range
                #Iterate over forward primers for this sequence
                current_fwd_primer = sequence.sequence_raw[forward_end_index - primer_width - offset : forward_end_index - offset]
                if calculate_degeneracy(current_fwd_primer) <= 4**5: #only proceed if the sequence is not "too degenerate"
                    for forward_primer in disambiguate(current_fwd_primer):
                        if forward_primer in self.primer2index['forward']:
                            if self.index2primer['forward'][self.primer2index['forward'][forward_primer]].feasible:
                                amplicon.primers['forward'][sequence.id].add(self.primer2index['forward'][forward_primer])
                                amplicon.full_primerset['forward'].add(self.primer2index['forward'][forward_primer])
                #Iterate over reverse primers for this sequence        
                current_rev_primer = reverse_complement(sequence.sequence_raw[reverse_start_index + offset : reverse_start_index + primer_width + offset])
                if calculate_degeneracy(current_rev_primer) <= 4**5:
                    for reverse_primer in disambiguate(current_rev_primer):
                        if reverse_primer in self.primer2index['reverse']:
                            if self.index2primer['reverse'][self.primer2index['reverse'][reverse_primer]].feasible:
                                amplicon.primers['reverse'][sequence.id].add(self.primer2index['reverse'][reverse_primer])
                                amplicon.full_primerset['reverse'].add(self.primer2index['reverse'][reverse_primer])

    def update_conflict_matrix(self, primers):
        '''
        Function that generates the conflicts for all the primer pairs that can be obtained by taking combinations of primers from $primers. If this PrimerIndex
        already has a conflict matrix, it will only be updated and not generated again.

        Parameters
        ----------
        primers : list[ Primer ]
            List of Primer objects to determine conflicts between.

        Returns
        -------
        None.

        '''
        if not self.conflict_matrix:
            #Matrix entry will be equal to -1 if not yet assigned, 1 if primers have a conflict, 2 if primers don't have a conflict
            self.conflict_matrix = {('f','r') : -1 * np.ones((len(self.index2primer['forward']), len(self.index2primer['reverse'])), dtype=np.int8), 
                                    ('f','f') : -1 * np.ones((len(self.index2primer['forward']), len(self.index2primer['forward'])), dtype=np.int8),
                                    ('r','r') : -1 * np.ones((len(self.index2primer['reverse']), len(self.index2primer['reverse'])), dtype=np.int8)}
        #Iterate over primer pairs
        for pair in itertools.combinations(primers, 2):
            #First primer is forward, second is reverse
            if (pair[0].orientation == 'forward' and pair[1].orientation == 'reverse'):
                current_index_pair = (self.primer2index['forward'][pair[0].sequence], self.primer2index['reverse'][pair[1].sequence])
                if self.conflict_matrix[('f','r')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix, self.thresholds['self_complementarity_threshold'])[0] > self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('f','r')][current_index_pair] = 1
                    else:
                        self.conflict_matrix[('f','r')][current_index_pair] = 2
            #First primer is reverse, second is forward
            elif (pair[0].orientation == 'reverse' and pair[1].orientation == 'forward'):
                current_index_pair = (self.primer2index['forward'][pair[1].sequence], self.primer2index['reverse'][pair[0].sequence])
                if self.conflict_matrix[('f','r')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix, self.thresholds['self_complementarity_threshold'])[0] > self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('f','r')][current_index_pair] = 1
                    else:
                        self.conflict_matrix[('f','r')][current_index_pair] = 2
            #Both primers are forward
            elif (pair[0].orientation == 'forward' and pair[1].orientation == 'forward'):
                current_index_pair = (self.primer2index['forward'][pair[0].sequence], self.primer2index['forward'][pair[1].sequence])
                if self.conflict_matrix[('f','f')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix, self.thresholds['self_complementarity_threshold'])[0] > self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('f','f')][current_index_pair] = 1
                        self.conflict_matrix[('f','f')][current_index_pair[1],current_index_pair[0]] = 1
                    else:
                        self.conflict_matrix[('f','f')][current_index_pair] = 2
                        self.conflict_matrix[('f','f')][current_index_pair[1],current_index_pair[0]] = 2
            #Both primers are reverse
            else:
                current_index_pair = (self.primer2index['reverse'][pair[0].sequence], self.primer2index['reverse'][pair[1].sequence])
                if self.conflict_matrix[('r','r')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix, self.thresholds['self_complementarity_threshold'])[0] > self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('r','r')][current_index_pair] = 1
                        self.conflict_matrix[('r','r')][current_index_pair[1],current_index_pair[0]] = 1
                    else:
                        self.conflict_matrix[('r','r')][current_index_pair] = 2
                        self.conflict_matrix[('r','r')][current_index_pair[1],current_index_pair[0]] = 2
                        
    def check_conflict(self, primer_pair):
        self.update_conflict_matrix(primer_pair)
        if primer_pair[0].orientation == 'forward' and primer_pair[1].orientation == 'reverse':
            orientation = ('f','r')
            pair = (self.primer2index['forward'][primer_pair[0].sequence], self.primer2index['reverse'][primer_pair[1].sequence])
        elif primer_pair[1].orientation == 'reverse' and primer_pair[1].orientation == 'forward':
            orientation = ('f','r')
            pair = (self.primer2index['forward'][primer_pair[1].sequence], self.primer2index['reverse'][primer_pair[0].sequence])
        else:
            orientation = (primer_pair[0].orientation[0], primer_pair[1].orientation[0])
            pair = (self.primer2index[primer_pair[0].orientation][primer_pair[0].sequence], self.primer2index[primer_pair[1].orientation][primer_pair[1].sequence])
        return self.conflict_matrix[orientation][pair]

    @staticmethod
    def generate_index(sequences, width, comparison_matrix):
        '''
        Static function that generates a primer index for the given sequences using a primer width of $width. For the
        multiprocessing variant see generate_index_mp

        Parameters
        ----------
        sequences : list[ Sequence ]
            List of sequences to find primers in.
        width : int
            Width of the primers to include in this index.
        comparison_matrix : dict[ (char,char) ]
            Dictionary that determines which characters should be considered equal.

        Returns
        -------
        primer_index : PrimerIndex
            Primer index containing all the primers of length #width that appear in the sequences in $sequences.

        '''
        primer_index = PrimerIndex()
        i = 0
        st = time.time()
        if type(sequences) == list: #If multiple sequences supplied
            for sequence in sequences:
                for cur_index in range(sequence.length_raw - width + 1):
                    current_fwd_primer = sequence.sequence_raw[cur_index : cur_index + width]
                    if calculate_degeneracy(current_fwd_primer) <= 4**5:
                        for forward_primer in disambiguate(current_fwd_primer):
                            primer_index.add_sequence(sequence, cur_index, forward_primer, 'forward')
                    current_rev_primer = reverse_complement(current_fwd_primer)
                    if calculate_degeneracy(current_rev_primer) <= 4**5:
                        for reverse_primer in disambiguate(current_rev_primer):
                            primer_index.add_sequence(sequence, cur_index, reverse_primer, 'reverse')
                i += 1
        else: #If a singular sequence is supplied
            for cur_index in range(sequences.length_raw - width + 1):
                current_fwd_primer = sequences.sequence_raw[cur_index : cur_index + width]
                if calculate_degeneracy(current_fwd_primer) <= 4**5:
                    for forward_primer in disambiguate(current_fwd_primer):
                        primer_index.add_sequence(sequences, cur_index, forward_primer, 'forward')
                current_rev_primer = reverse_complement(current_fwd_primer)
                if calculate_degeneracy(current_rev_primer) <= 4**5:
                    for reverse_primer in disambiguate(current_rev_primer):
                        primer_index.add_sequence(sequences, cur_index, reverse_primer, 'reverse')
            i += 1
        print('Time elapsed processing %d sequences: %fs' % (i,(time.time() - st)))
        return primer_index
    
    @staticmethod
    def generate_index_mp(sequences, width, comparison_matrix, processors=1):
        if processors > 1:
            sequences_partitioned = [ sequences[i:i+(ceil(len(sequences)/processors))] for i in range(0, len(sequences), ceil(len(sequences)/processors))]
            with mp.Pool(processors) as pool:
                indices = pool.starmap(PrimerIndex.generate_index, zip(sequences_partitioned, itertools.repeat(width), itertools.repeat(comparison_matrix)))
            master_index = indices[0]
            for index in indices[1:]:
                master_index.merge_indices(index)
            return master_index
        else:
            return PrimerIndex.generate_index(sequences, width, comparison_matrix)
        
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
    args = parser.parse_args()
    
    
    #Initialize variables to store information
    runtimes = []
    cur_time = time.time()
    
    #Read sequences
    st = time.time()
    sequences = generate_sequences(args.metadata, args.sequences)
    with open(args.output + '/runtimes.txt', 'a') as f:
        f.write('Time spent generating sequences: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent generating sequences: ' + str(time.time() - st))
    
    #Generate comparison matrix
    comparison = generate_opportunistic_matrix()
    
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
        sequences, lb, ub, feasible_amplicons, _ = preprocess_sequences(sequences, args.search_width, variants_location=args.variants_location, variants=variants, amplicon_width=args.amplicon_width, misalign_threshold=args.misalign_threshold)
        with open(args.output + '/runtimes.txt', 'a') as f:
            f.write('Variants considered:\n')
            for variant in variants:
                f.write(variant + '\n')
            f.write('Total sequences = ' + str(len(sequences)) + '\n')
    else:
        sequences, lb, ub, feasible_amplicons, _ = preprocess_sequences(sequences, args.search_width, amplicon_width=args.amplicon_width, misalign_threshold=args.misalign_threshold)
        with open(args.output + '/runtimes.txt', 'a') as f:
            f.write('Variants considered:\n')
            for variant in ['Alpha','Beta','Gamma','Delta','Epsilon','Zeta','Eta','Kappa','Mu','Omicron']:
                f.write(variant + '\n')
            f.write('Total sequences = ' + str(len(sequences)) + '\n')
    
    with open(args.output + '/runtimes.txt', 'a') as f:
        f.write('Time spent pre-processing sequences and determining feasible amplicons: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent pre-processing sequences and determining feasible amplicons: ' + str(time.time() - st))
    
    #Generate primer index
    st = time.time()
    PI = PrimerIndex.generate_index_mp(sequences, args.primer_width, comparison, processors=args.cores)
    PI.remove_redundant()
    with open(args.output + '/runtimes.txt', 'a') as f:
        f.write('Time spent generating primer index and filtering for feasible primers: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent generating primer index and filtering for feasible primers: ' + str(time.time() -st))
    
    #Generate amplicons
    st = time.time()
    amplicons = generate_amplicons_mp_exp(sequences, args.amplicon_width, comparison, feasible_amplicons=feasible_amplicons, processors=args.cores, amplicon_threshold=args.amplicon_threshold)
    with open(args.output + '/runtimes.txt', 'a') as f:
        f.write('Time spent generating amplicon differentiation ' + str(time.time() - st) + '\n')
        f.write('Total feasible amplicons: ' + str(len(amplicons)) + '\n')
    #runtimes.append('Time spent generating amplicon differentiation: ' + str(time.time() - st))
    
    #Run greedy
    st = time.time()
    logs, amplicons, result_amplicons = greedy(sequences, amplicons, args.primer_width, args.search_width, PI, comparison, args.amplicons, args.coverage, 5, logging=True, multiplex=args.multiplex)
    with open(args.output + '/runtimes.txt', 'a') as f:
        f.write('Time spent running greedy algorithm: ' + str(time.time() - st) + '\n')
    #runtimes.append('Time spent running greedy algorithm: ' + str(time.time() - st))
    
    #Run final optimization
    if args.multiplex:
        st = time.time()
        check_primer_feasibility(sequences, result_amplicons, PI, optimize=1, coverage=args.coverage)
        with open(args.output + '/runtimes.txt', 'a') as f:
            f.write('Time spent doing final primer optimization: ' + str(time.time() - st))
    else:
        st = time.time()
        for amplicon in result_amplicons:
            cur_primers = check_primer_feasibility(sequences, [amplicon], PI, optimize=1, coverage=args.coverage)
            with open(args.output + '/runtimes.txt', 'a') as f:
                f.write('Amplicon: ' + str(amplicon.id) + '\n')
                f.write('Forward primers\n')
                for fwd in cur_primers['forward']:
                    f.write(fwd + '\n')
                f.write('Reverse primers\n')
                for rev in cur_primers['reverse']:
                    f.write(rev + '\n')
                
    #runtimes.append('Time spent doing final primer optimization: ' + str(time.time() - st))
    
    #with open(args.output + '/runtimes.txt', 'w') as f:
    #    for line in runtimes:
    #        f.write(line + '\n')
    with open(args.output + '/logfile.txt', 'w') as f:
        for line in logs:
            f.write(line + '\n')
    
    # primer_width = 25
    # search_width = 100
    # amplicon_width = 200
    # misalign_threshold = 2
    
    # sequences = generate_sequences('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/data/Global/global_all_time_N0','/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/data/Global/global_all_time_N0/unfiltered')
    # st = time.time()
    # sequences, lb, ub, feasible_amplicons, _ = preprocess_sequences(sequences, search_width, amplicon_width=amplicon_width, misalign_threshold=misalign_threshold)
    # print( 'Time determining feasible amplicons and preprocessing sequences: %fs' % (time.time() - st) )
    
    # st = time.time()
    # PI = PrimerIndex.generate_index_mp(sequences[:500], 25, generate_opportunistic_matrix(), processors=8)
    # PI.remove_redundant()
    # print( 'Time spent generating new index using MP: %fs' % (time.time() - st) )
    
    # st = time.time()
    # amplicons = list(feasible_amplicons)
    # amplicons.sort()
    # A = generate_amplicons_mp(sequences[:500], amplicon_width, generate_opportunistic_matrix(), feasible_amplicons=amplicons, processors=8, amplicon_threshold=0)
    # print( 'Time spent generating amplicons using MP: %fs' % (time.time() - st) )
    
    # A = Amplicon(4000, 4200)
    # st = time.time()
    # primers_unfiltered = PI_mp.check_amplicon(sequences, A, 25, 50)
    # print( 'Time spent finding primers in unfiltered index: %fs' % (time.time() - st) )

    # st = time.time()
    # PI_mp.remove_redundant()
    # print( 'Time spent filtering primer index: %fs' % (time.time() - st) )
    
    # st = time.time()
    # primers_filtered = PI_mp.check_amplicon(sequences, A, 25, 50)
    # print( 'Time spent finding primers in unfiltered index: %fs' % (time.time() - st) )
            
"""
class Sequence:
    def __init__(self, sequence, identifier, lineage=None):
        self.sequence = sequence
        self.length = len(sequence)
        self.sequence_raw = sequence.replace('-','')
        self.length_raw = len(self.sequence_raw)
        self.id = identifier
        self.lineage = lineage #can be None
        self.aligned_to_trim = [] #Initially empty
        self.kmers = {'forward' : {}, 'reverse' : {}}
        
    def __eq__(self, other):
        try:
            return self.id == other.id
        except:
            return False
    
    def __hash__(self):
        return hash(self.id)
    
    def __repr__(self):
        return 'Sequence'
    
    def reverseComplement(self, rev=True):
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
        for i in range(len(self.sequence_raw)):
            res += translate[self.sequence_raw[i]] #complement every nucleotide
        if rev:
            return res[::-1]
        else:
            return res
        
    def alignToTrim(self):
        self.aligned_to_trim = [0]*self.length
        map_index = -1
        for char_index in range(self.length):
            if self.sequence[char_index] != '-':
                map_index += 1
            self.aligned_to_trim[char_index] = map_index
        return self.aligned_to_trim
    
    def findBounds(self, min_non_align):
        lb = self.aligned_to_trim.index(min_non_align)
        ub = self.aligned_to_trim.index(self.length_raw - min_non_align)
        return lb, ub
        
    def generateKmers(self, k, comparison_matrix):
        '''
        Function that generates all the kmers of this sequence, and their location in the sequence

        Parameters
        ----------
        k : int
            Value for the k in k-mer.
        comparison_matrix : dict[(char,char)]
            Dictionary that contains information about character comparisons.

        Returns
        -------
        self.kmers
            k-mer dictionary created in the process.

        '''
        self.kmers = {'forward' : {}, 'reverse' : {}}
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
        
        reverse_sequence = self.reverseComplement()
        for start_index in range(len(self.sequence_raw) - k + 1):
            current_kmer_fwd = self.sequence_raw[start_index:start_index+k]
            current_kmer_rev = reverse_sequence[start_index:start_index+k]
            if calcDegeneracy(current_kmer_fwd) <= 4**5: #Check if current k-mar is not "too degenerate"
                current_disambiguated = disambiguate(current_kmer_fwd)
                for kmer_candidate in current_disambiguated:
                    if kmer_candidate in self.kmers['forward']:
                        self.kmers['forward'][kmer_candidate].append((start_index, start_index+k))
                    else:
                        self.kmers['forward'][kmer_candidate] = [(start_index, start_index+k)]
            if calcDegeneracy(current_kmer_rev) <= 4**5:
                current_disambiguated = disambiguate(current_kmer_rev)
                for kmer_candidate in current_disambiguated:
                    if kmer_candidate in self.kmers['reverse']:
                        self.kmers['reverse'][kmer_candidate].append((len(self.sequence_raw) - start_index - k, len(self.sequence_raw) - start_index))
                    else:
                        self.kmers['reverse'][kmer_candidate] = [(len(self.sequence_raw) - start_index - k, len(self.sequence_raw) - start_index)]
        return self.kmers   
    
    def generateKmers_exp(self, k, comparison_matrix, illegal_kmers=[]):
        '''
        Function that generates all the kmers of this sequence, and their location in the sequence

        Parameters
        ----------
        k : int
            Value for the k in k-mer.
        comparison_matrix : dict[(char,char)]
            Dictionary that contains information about character comparisons.

        Returns
        -------
        self.kmers
            k-mer dictionary created in the process.

        '''
        self.kmers = {'forward' : {}, 'reverse' : {}}
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
        reverse_sequence = self.reverseComplement()
        for start_index in range(len(self.sequence_raw) - k + 1):
            current_kmer_fwd = self.sequence_raw[start_index:start_index+k]
            current_kmer_rev = reverse_sequence[start_index:start_index+k]
            if calcDegeneracy(current_kmer_fwd) <= 4**5: #Check if current k-mer is not "too degenerate"
                current_disambiguated = disambiguate(current_kmer_fwd)
                for kmer_candidate in current_disambiguated:
                    if kmer_candidate not in illegal_kmers:
                        tmp = Primer(kmer_candidate, None, 'forward')
                        if tmp.determineProperties(comparison_matrix, 1000, 1000):
                            if kmer_candidate in self.kmers['forward']:
                                self.kmers['forward'][kmer_candidate].append((start_index, start_index+k))
                            else:
                                self.kmers['forward'][kmer_candidate] = [(start_index, start_index+k)]
                        else:
                            illegal_kmers.append(kmer_candidate)
            if calcDegeneracy(current_kmer_rev) <= 4**5:
                current_disambiguated = disambiguate(current_kmer_rev)
                for kmer_candidate in current_disambiguated:
                    if kmer_candidate not in illegal_kmers:
                        tmp = Primer(kmer_candidate, None, 'reverse')
                        if tmp.determineProperties(comparison_matrix, 1000, 1000):
                            if kmer_candidate in self.kmers['reverse']:
                                self.kmers['reverse'][kmer_candidate].append((len(self.sequence_raw) - start_index - k, len(self.sequence_raw) - start_index))
                            else:
                                self.kmers['reverse'][kmer_candidate] = [(len(self.sequence_raw) - start_index - k, len(self.sequence_raw) - start_index)]
                        else:
                            illegal_kmers.append(kmer_candidate)
        return illegal_kmers 
    
    def toText(self):
        res = self.id + '\n' + self.lineage + '\n' + self.sequence + '\n'
        for index in self.aligned_to_trim:
            res += str(index) + ','
        res = res[:-1]
        res += '\n' + str(self.kmers) + '\n'
        return res
    
    @staticmethod
    def readSequences(file, n_cores=None):
        def processSequence(in_lines):
            res = Sequence(in_lines[2].rstrip(), in_lines[0].rstrip(), lineage=in_lines[1].rstrip())
            res.aligned_to_trim = list(eval(in_lines[3].rstrip()))
            res.kmers = eval(in_lines[4].rstrip())
            return res
        res = []
        with open(file, 'r') as f:
            lines = f.readlines()
        if n_cores:
            actual_lines = []
            for i in range(0, len(lines), 5):
                actual_lines.append([lines[i+k] for k in range(5)])
            pool = mp.Pool(n_cores)
            res = pool.map(processSequence, actual_lines)
            pool.close()
        else:
            for i in range(0, len(lines), 5):
                res.append(Sequence(lines[i+2].rstrip(), lines[i].rstrip(), lineage=lines[i+1].rstrip()))
                res[-1].aligned_to_trim = list(eval(lines[i+3].rstrip()))
                res[-1].kmers = eval(lines[i+4].rstrip())
        return res
        
class Amplicon:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.id = (start,end)
        self.primers = {'forward' : [], 'reverse' : []}
        self.primers_per_sequence = {}
        self.differences = set()
        self.differences_proper = set()
        self.weight = 0
        
    def __eq__(self, other):
        try:
            return self.id == other.id
        except:
            return False
    
    def __lt__(self, other):
        if type(other) == Amplicon:
            return self.start < other.start
        
    def __repr__(self):
        return 'Amplicon'
    
    def addPrimers(self, primers, orientation):
        if type(primers) == list:
            for primer in primers:
                self.primers[primer.orientation].append(primer)
                for sequence in primer.binding_sequences:
                    if sequence in self.primers_per_sequence:
                        self.primers_per_sequence[sequence][primer.orientation].append(primer)
                    else:
                        self.primers_per_sequence[sequence] = {'forward' : [], 'reverse' : []}
                        self.primers_per_sequence[sequence][primer.orientation].append(primer)
        else:
            self.primers[primers.orientation].append(primers)
            for sequence in primers.binding_sequences:
                if sequence in self.primers_per_sequence:
                    self.primers_per_sequence[sequence][primers.orientation].append(primers)
                else:
                    self.primers_per_sequence[sequence] = {'forward' : [], 'reverse' : []}
                    self.primers_per_sequence[sequence][primers.orientation].append(primers)
            
    def setDifferences(self, differences):
        self.differences = set(differences)
        for difference in differences:
            if difference[0].lineage != difference[1].lineage:
                self.differences_proper.add(difference)
        
    def setWeight(self, lineages, ignore, weighted=False):
        def calcWeight(lineage1, lineage2):
            lin1_split = lineage1.split('.')
            lin2_split = lineage2.split('.')
            for i in range(len(lin1_split)):
                if i >= len(lin2_split):
                    return len(lin1_split) - len(lin2_split)
                elif lin1_split[i] != lin2_split[i]:
                    return len(lin1_split[i:]) + len(lin2_split[i:])
            return len(lin2_split) - len(lin1_split)
        
        self.weight = 0
        if weighted:
            for pair in (self.differences - ignore):
                self.weight += calcWeight(lineages[pair[0]], lineages[pair[1]])
        else:
            self.weight = len(self.differences - ignore)
        return self.weight
    
    def toText(self):
        res = str(self.start) + ',' + str(self.end) + '\n'
        
        if len(self.differences) > 0:
            res += '{'
            for difference in self.differences:
                res += "('" + difference[0].id + "', '" + difference[1].id + "'), "
            res = res[:-2] + '}\n'
        else:
            res + '{}\n'
        '''
        if len(self.differences_proper) > 0:
            res += '{'
            for difference in self.differences_proper:
                res += "('" + difference[0].id + "', '" + difference[1].id + "'), "
            res = res[:-2] + '}\n'
        else:
            res + '{}\n'
        ''' 
        res += "{'forward': ["
        if len(self.primers['forward']) > 0:
            for primer in self.primers['forward']:
                res += "'" + primer.sequence + "', "
            res = res[:-2] + '], '
        else:
            res += '], '
        res += "'reverse': ["
        if len(self.primers['reverse']) > 0:
            for primer in self.primers['reverse']:
                res += "'" + primer.sequence + "', "
            res = res[:-2] + ']}\n'
        else:
            res += ']}\n'
            
        res += '{'
        if len(self.primers_per_sequence) > 0:
            for sequence in self.primers_per_sequence:
                res += "'" + sequence.id + "': {'forward': ["
                if len(self.primers_per_sequence[sequence]['forward']) > 0:
                    for primer in self.primers_per_sequence[sequence]['forward']:
                        res += "'" + primer.sequence + "', "
                    res = res[:-2] + "], 'reverse': ["
                else:
                    res += "], 'reverse: ["
                if len(self.primers_per_sequence[sequence]['reverse']) > 0:
                    for primer in self.primers_per_sequence[sequence]['reverse']:
                        res += "'" + primer.sequence + "', "
                    res = res[:-2] + "]}, "
                else:
                    res += ']}, '
            res = res[:-2] + '}\n'
        else:
            res += '}\n'
        return res
    
    @staticmethod
    def readAmplicons(file, sequences, primers):
        def processAmplicon(in_lines, sequences, primers):
            res = Amplicon(int(in_lines[0].rstrip().split(',')[0]), int(in_lines[0].rstrip().split(',')[1]))
            for difference in eval(in_lines[1]):
                proxy_sequence_1 = Sequence('n', difference[0])
                proxy_sequence_2 = Sequence('n', difference[1])
                cur_diff = (sequences[sequences.index(proxy_sequence_1)], sequences[sequences.index(proxy_sequence_2)])
                res.differences.add(cur_diff)
                if cur_diff[0].lineage != cur_diff[1].lineage:
                    res.differences_proper.add(cur_diff)
            return res
        
        with open(file, 'r') as f:
            lines = f.readlines()
        res = []
        for i in range(0, len(lines), 4):
            res.append(processAmplicon(lines[i:i+4], sequences, primers))
        return res
            
        

class Primer:
    def __init__(self, sequence, amplicon, orientation, comparison_tolerance=13):
        self.sequence = sequence
        self.amplicon = amplicon
        self.orientation = orientation
        self.binding_sequences = []
        self.comparison_tolerance = comparison_tolerance
        
    def __eq__(self, other):
        try:
            return self.sequence == other.sequence
        except:
            return False
    
    def __hash__(self):
        return hash(self.sequence)
    
    def __repr__(self):
        return 'Primer'
    
    def reverseComplement(self, rev=True):
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
        for i in range(len(self.sequence)):
            res += translate[self.sequence[i]] #complement every nucleotide
        if rev:
            return res[::-1]
        else:
            return res
    
    def checkCompatibility(self, other, comparison_matrix):
        res = (0, '')
        
        to_compare = other.reverseComplement()
        
        for i in range(len(self.sequence)):
            fwd_check = []
            rev_check = []
            for j in range(i, len(self.sequence)):
                #Check forward comparison
                if comparison_matrix[(to_compare[j], self.sequence[j-i])][0] and not comparison_matrix[(to_compare[j], self.sequence[j-i])][1]:
                    fwd_check.append(1)
                else:
                    fwd_check.append(0)
                #Check backward comparison
                if comparison_matrix[(self.sequence[j], to_compare[j-i])][0] and not comparison_matrix[(self.sequence[j], to_compare[j-i])][1]:
                    rev_check.append(1)
                else:
                    rev_check.append(0)
            comparison = sum(fwd_check)
            if comparison > res[0]:
                if self == other:
                    res = [comparison, i*'*' + self.sequence + '\n' + self.sequence[::-1] + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in fwd_check) + i*'*']
                else:
                    res = [comparison, i*'*' + self.sequence + '\n' + other.sequence[::-1] + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in fwd_check) + i*'*']
            comparison = sum(rev_check)
            if comparison > res[0]:
                if self == other:
                    res = [comparison, i*'*' + self.sequence[::-1] + '\n' + self.sequence + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in rev_check) + i*'*']
                else:
                    res = [comparison, i*'*' + other.sequence[::-1] + '\n' + self.sequence + i*'*' + '\n' + i*'*' + ''.join(str(c) for c in rev_check) + i*'*']
            #If number of matches exceeds upperbound (10 so < 11) then return
            if comparison >= self.comparison_tolerance:
                return res
        return res

    def determineProperties(self, comparison_matrix, at_threshold, gc_threshold, log=False):
        def calcDegeneracy(self):
            res = 1
            for char in self.sequence:
                if char in ['y','r','s','w','m','k']:
                    res = res*2
                elif char in ['b','d','h','v']:
                    res = res*3
                elif char == 'n':
                    res = res*4
            return res
        
        def calcGC(self):
            res = len(self.sequence)
            for char in self.sequence:
                if char in ['a','t','w','-']:
                    res -= 1
            return res/len(self.sequence)
        
        def calcEndStats(self, comparison_matrix, at_threshold, gc_threshold):
            res = [0, 0, False]
            for i in range(1,6):
                if i <= 3 and (comparison_matrix[(self.sequence[-i], 'a')][0] or comparison_matrix[(self.sequence[-i], 't')][0]):
                    res[0] += 1
                if comparison_matrix[(self.sequence[-i], 'c')][0] or comparison_matrix[(self.sequence[-i], 'g')][0]:
                    res[1] += 1
                    if i == 1:
                        res[2] = True
            return res
        
        def calcLongestMonoRun(self, comparison_matrix):
            longest = 0
            current = 1
            left_char = self.sequence[0]
            for right_char in self.sequence[1:]:
                if comparison_matrix[(left_char, right_char)][0]:
                    current += 1
                else:
                    longest = max(longest, current)
                    current = 1
                    left_char = right_char
            return max(longest, current)
        
        def calcLongestDuoRun(self, comparison_matrix):
            longest = 0
            current = 1
            left_char = self.sequence[0]
            right_char = self.sequence[1]
            index = 2
            while index < len(self.sequence) - 1:
                if comparison_matrix[(left_char, self.sequence[index])][0] and comparison_matrix[(right_char, self.sequence[index+1])][0]:
                    current += 1
                    index += 2
                else:
                    longest = max(longest, current)
                    current = 1
                    left_char = self.sequence[index-1]
                    right_char = self.sequence[index]
                    index += 1
            return max(longest, current)
        #Assign properties
        self.degeneracy = calcDegeneracy(self)
        self.gc_content = calcGC(self)
        try:
            self.melting_temperature = mt.Tm_NN(self.sequence, strict=False)
        except:
            self.melting_temperature = 0
        stats = calcEndStats(self, comparison_matrix, at_threshold, gc_threshold)
        self.gc_clamp = stats[2]
        self.end_at_count = stats[0]
        self.end_gc_count = stats[1]
        self.max_run = calcLongestMonoRun(self, comparison_matrix)
        self.max_2_run = calcLongestDuoRun(self, comparison_matrix)
        self.MFE = RNA.fold(self.sequence)[1]
        self.self_dimer = self.checkCompatibility(self, comparison_matrix)
        #Check if all properties are met for primers
        if self.degeneracy > 1:
            if log:
                print('Too degenerate')
            return False
        if self.sequence.count('-') > 0:
            if log:
                print('Too many misaligns')
            return False
        if self.gc_content < 0.4 or self.gc_content > 0.6:
            if log:
                print('GC-content')
            return False
        if self.melting_temperature < 55 or self.melting_temperature > 75:
            if log:
                print('melting temperature')
            return False
        if self.end_at_count >= at_threshold:
            if log:
                print('at content at 3-end too high')
            return False
        if self.end_gc_count >= gc_threshold:
            if log:
                print('gc content at 3-end too high')
            return False
        if self.max_run >= 4:
            if log:
                print('Too long monorun')
            return False
        if self.max_2_run >= 3:
            if log:
                print('Too long duorun')
            return False
        if self.MFE <= -5:
            if log:
                print('MFE too low')
            return False
        if self.self_dimer[0] >= self.comparison_tolerance: #If self-complementarity is over 10 nt-pairs then reject
            if log:
                print('Self dimer problem')
            return False
        return True

    def addBindingSequence(self, sequence):
        self.binding_sequences.append(sequence)
        
    def checkBindingEvents(self, sequences):
        for sequence in sequences:
            if self.sequence in sequence.kmers[self.orientation] and len(sequence.kmers[self.orientation][self.sequence]) == 1:
                self.binding_sequences.append(sequence)
        return set(sequences) == set(self.binding_sequences)
    
    def toText(self):
        res = self.sequence + '\n'
        res += str(self.amplicon.start) + ',' + str(self.amplicon.end) + '\n'
        res += self.orientation +'\n'
        res += '['
        for sequence in self.binding_sequences:
            res += "'" + sequence.id + "', "
        res = res[:-2] + ']\n'
        return res
    
    @staticmethod
    def readPrimers(file, sequences):
        def processPrimer(in_lines, sequences):
            res = Primer(in_lines[0], (int(in_lines[1].split(',')[0]), int(in_lines[1].split(',')[1])), in_lines[2])
            return res
        with open(file, 'r') as f:
            lines = f.readlines()
        res = []
        for i in range(0, len(lines), 4):
            res.append(processPrimer(lines[i:i+4], sequences))
        return res
"""      