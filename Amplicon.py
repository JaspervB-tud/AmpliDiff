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
                    self.differences_proper.add((difference[0].alt_id, difference[1].alt_id))
        except:
            pass
        
    def check_differences(self, sequences):
        for difference in self.differences:
            try:
                if sequences[sequences.index(difference[0])].lineage != sequences[sequences.index(difference[1])].lineage:
                    self.differences_proper.add((difference[0], difference[1]))
            except:
                continue