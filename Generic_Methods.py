import itertools

def equivalent_characters(c):
    if c == '-':
        return set(['-'])
    else:
        res = set()
        if c in ['a','r','m','w','d','h','v','n']:
            res.add('a')
        if c in ['c','y','m','s','b','h','v','n']:
            res.add('c')
        if c in ['g','r','k','s','b','d','v','n']:
            res.add('g')
        if c in ['t','y','k','w','b','d','h','n']:
            res.add('t')
        return res

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

def generate_comparison_table():
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
            if c1 == c2:
                res[(c1,c2)] = True
            elif c1 == 'n' and c2 != '-':
                res[(c1,c2)] = True
            elif c2 == 'n' and c2 != '-':
                res[(c1,c2)] = True
            elif c1 in char_comp[c2] and c2 in char_comp[c1]:
                res[(c1,c2)] = True
            else:
                res[(c1,c2)] = False
    return res

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