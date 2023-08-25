import itertools

########################### Standalone functions ###########################
def equivalent_characters(c):
    '''
    Function that returns a set of non-degenerate nucleotides that are equivalent to the (possibly degenerate) input nucleotide.

    Parameters
    ----------
    c : char
        Nucleotide according to IUPAC notation for which equivalent "base" nucleotides will be returned.

    Returns
    -------
    equiv : set
        Set of characters equivalent to input character according to IUPAC notation.

    '''
    if c == '-':
        return set(['-'])
    else:
        equiv = set()
        if c in ['a','r','m','w','d','h','v','n']:
            equiv.add('a')
        if c in ['c','y','m','s','b','h','v','n']:
            equiv.add('c')
        if c in ['g','r','k','s','b','d','v','n']:
            equiv.add('g')
        if c in ['t','y','k','w','b','d','h','n']:
            equiv.add('t')
        return equiv

def generate_comparison_matrix():
    '''
    Function that constructs a dictionary with pairs of IUPAC nucleotide characters as keys, and a tuple of
    with booleans (first element is True if characters share at least one non-degenerate nucleotide, and
    second element is True if either character is a misalignment character) as keys.

    Returns
    -------
    res : dict[ (char,char) ] -> (bool, bool)
        Dictionary with character tuples as keys, and boolean tuples as values.

    '''
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
    comparison_matrix = {}
    for c1 in chars:
        for c2 in chars:
            if c1 == '-' or c2 == '-':
                if c1 == c2:
                    comparison_matrix[(c1,c2)] = (True,True)
                else:
                    comparison_matrix[(c1,c2)] = (False,True)
            elif c1 == 'n' or c2 == 'n':
                comparison_matrix[(c1,c2)] = (True,False)
            elif c1 in char_comp[c2] and c2 in char_comp[c1]:
                comparison_matrix[(c1,c2)] = (True,False)
            else:
                comparison_matrix[(c1,c2)] = (False,False)
    return comparison_matrix

def reverse_complement(sequence, rev=True):
    '''
    Function that returns the (reverse) complement of the given sequence.

    Parameters
    ----------
    sequence : str
        String representation of the sequence to determine reverse complement of.
    rev : bool, optional
        True if the reverse should be returned, False otherwise. The default is True.

    Returns
    -------
    rev_comp : str
        (Reverse) complement of the input sequence.

    '''
    #Define the complement of every possible IUPAC nucleotide
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
    rev_comp = ''
    for i in range(len(sequence)):
        rev_comp += translate[sequence[i]]
    if rev:
        return rev_comp[::-1]
    else:
        return rev_comp

def disambiguate(sequence):
    '''
    Function that disambiguates the given sequence by generating all its non-degenerate representations.

    Parameters
    ----------
    sequence : str
        String representation of the sequence to disambiguate.

    Returns
    -------
    repr : [str]
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
    
    repr = translation[sequence[0]].copy()
    for char in sequence[1:]:
        for subsequence_index in range(len(repr)):
            new_subsequences = []
            for new_char in translation[char]:
                new_subsequences.append(repr[subsequence_index] + new_char)
            repr[subsequence_index] = new_subsequences
        repr = list(itertools.chain(*repr))
    return repr

########################### Calculation functions for sequence properties ###########################

def calculate_degeneracy(sequence):
    '''
    Function that returns the degeneracy of a sequence of nucleotides which is defined as the
    cardinality of the set of all possible non-degenerate representations of the sequence.

    Parameters
    ----------
    sequence : str
        String representation of a series of consecutive nucleotides.

    Returns
    -------
    degen : int
        Degeneracy of the input sequence.

    '''
    degen = 1
    for char in sequence:
        if char in ['y', 'r', 's', 'w', 'm', 'k']:
            degen = degen*2
        elif char in ['b', 'd', 'h', 'v']:
            degen = degen*3
        elif char == 'n':
            degen = degen*4
    return degen

def calculate_GC(sequence):
    '''
    Function that calculates the GC-content of a (possibly degenerate) sequence.

    Parameters
    ----------
    sequence : str
        String representation of sequence of interest.

    Returns
    -------
    float
        GC-content of input sequence.

    '''
    gc = len(sequence)
    for char in sequence:
        if char in ['a','t','w','-']: #assumption is that any character that can represent a G or C counts here
            gc -= 1
    return gc / len(sequence)

def calculate_end_stats(sequence, comparison_matrix):
    '''
    Function that determines the number of a/t (c/g) characters in the last 3 (5) characters of the 3'-end of the input sequence.

    Parameters
    ----------
    sequence : str
        String representation of sequence of interest.
    comparison_matrix : dict[ (char,char) ]
        Dictionary that determines which characters should be considered equal.

    Returns
    -------
    stats : (int, int, bool)
        Triplet where the first element is the number of a/t chars in final 3, 
        second element the number of c/g in final 5 and last element is true when last character is c/g.

    '''
    stats = [0, 0, False]
    for i in range(1, 4):
        if comparison_matrix[(sequence[-i], 'a')][0] or comparison_matrix[(sequence[-i], 't')][0]:
            stats[0] += 1
        elif comparison_matrix[(sequence[-i], 'c')][0] or comparison_matrix[(sequence[-i], 'g')][0]:
            stats[1] += 1
            if i == 1:
                stats[2] = True
    for i in range(4, 6):
        if comparison_matrix[(sequence[-i], 'c')][0] or comparison_matrix[(sequence[-i], 'g')][0]:
            stats[1] += 1
    return stats

def calculate_longest_monorun(sequence, comparison_matrix):
    '''
    Function that calculates the longest run of a single character in the given sequence.

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
    Function that calculates the longest run of a pair of characters in the given sequence.

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