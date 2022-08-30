import numpy as np

class Sequence:
    existing_sequences = 0
    lineage_to_number = {}
    def __init__(self, sequence, identifier, lineage=None):
        self.sequence = sequence
        self.length = len(sequence)
        self.sequence_raw = sequence.replace('-','')
        self.length_raw = len(self.sequence_raw)
        self.id = identifier
        self.id_num = Sequence.existing_sequences
        Sequence.existing_sequences += 1
        self.lineage = lineage
        if lineage:
            if lineage in Sequence.lineage_to_number:
                self.lineage_num = Sequence.lineage_to_number[lineage]
            else:
                self.lineage_num = len(Sequence.lineage_to_number)
                Sequence.lineage_to_number[lineage] = len(Sequence.lineage_to_number)
        self.aligned_to_trim = np.zeros((1))
        
    def __eq__(self, other):
        try:
            if type(other) == str:
                return self.sequence == other or self.sequence_raw == other or self.id == other
            elif type(other) == int:
                return self.alt_id == other
            else:
                return self.id == other.id or self.alt_id == other.alt_id
        except:
            return False
        
    def __hash__(self):
        return hash(self.alt_id)
    
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