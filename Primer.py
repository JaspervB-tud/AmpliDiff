from classless_methods import reverse_complement, calculate_degeneracy, calculate_GC, calculate_end_stats, calculate_longest_monorun, calculate_longest_duorun
from Bio.SeqUtils import MeltingTemp as mt
import RNA

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
            if sequence.id_num not in self.indices:
                self.indices[sequence.id_num] = sequence_index
            else:
                if not self.indices[sequence.id_num] == sequence_index:
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
            Maximal number of allowed C/G characters in final 5 nucleotides (3'-end). The default is 3.
        monorun_threshold : int, optional
            Maximal length of a run of a single nucleotide. The default is 3.
        duorun_threshold : int, optional
            Maximal length of a run of a pair of nucleotides. The default is 2.
        mfe_threshold : float, optional
            Minimal allowed MFE representing the "risk" of hairpin formation. The default is -5.
        self_complementarity_threshold : int, optional
            Maximal number of complementary base-pairs in every alignment of this primer with its own reverse. The default is 10.
        verbose : bool, optional
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
            self.mfe_hairpin = RNA.fold(self.sequence)[1]
            if self.mfe_hairpin < mfe_threshold:
                problems.append('Risk of hairpin/secondary structure formation too high')
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
            if RNA.fold(self.sequence)[1] < mfe_threshold:
                self.feasible = False
                return False
            return True