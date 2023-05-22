from Bio import SeqIO
import numpy as np
from Sequence import *
from Generic_Methods import reverse_complement, calculate_degeneracy, generate_opportunistic_matrix
import math
import itertools
import csv
import argparse

def generate_sequences(seq_path, meta_path, max_n=10**10):
    '''
    Function that reads sequences from a fasta file and saves them as Sequence objects with metadata from "metadata.tsv".

    Parameters
    ----------
    seq_path: str
        Absolute path of the sequences.
    meta_path : str
        Absolute path to the metadata.tsv file
    max_n : int, optional
        Maximum number of fully degenerate nucleotides allowed. The default is 10.

    Returns
    -------
    sequences : list[Sequence]
        List of sequences contained in seq_path.

    '''
    sequences_temp = {}
    to_delete = []
    
    #Read sequences from fasta file
    sequences_object = SeqIO.parse(open(seq_path), 'fasta')
    for sequence in sequences_object:
        sequences_temp[sequence.id.split('|')[0]] = str(sequence.seq.lower())
        if len(sequence.seq.replace('-','')) < 29000 or sequence.seq.lower().count('n') > max_n: #require sequences to have length at least 29k
            to_delete.append(sequence.id.split('|')[0]) #store sequences that should be removed due to too many 'n's or being too short
    #Read metadata unless impossible in which case we assign every sequence its own "lineage"
    skip = -1
    num_processed = 0
    sequences = []
    for meta in csv.reader(open(meta_path), delimiter='\t'):
        try:
            if skip == -1: #first line is always the header line
                for cur_meta in range(len(meta)):
                    if 'lineage' in meta[cur_meta].lower():
                        skip = cur_meta
                        break
            else:
                #meta[0] contains id of the sequence
                if meta[0] not in to_delete:
                    sequences.append(Sequence(sequences_temp[meta[0].replace(' ','')], meta[0], lineage=meta[skip]))
                    num_processed += 1
        except: #unable to read metadata
            '''
            print('Unable to read metadata from file, making up lineages for every sequence')
            sequences = []
            i = 0
            for identifier in sequences_temp:
                sequences.append(Sequence(sequences_temp[identifier], identifier, lineage=str(i)))
                i += 1
            '''
            #print('id:', meta, 'not found, skipping')
            continue
    print('Number of sequences processed:', num_processed)
    return sequences

def read_bedfile(filename):
    '''
    Function that reads the primer.bed file of the ARTIC protocol.

    Parameters
    ----------
    filename : str
        Absolute path to the bed file.

    Returns
    -------
    amplicons : dict
        Dictionary with amplicon ids as keys and a dictionary {'forward':[],'reverse':[]} as value storing the primers per amplicon.
    primerlist : dict
        Dictionary with keys forward and reverse which are simply aggregates of the forward and reverse primers of every amplicon in $amplicons.

    '''
    amplicons = {}
    primerlist = {'forward': [], 'reverse': []}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.split('\t')
            amplicon_id = line[3].split('_')
            if amplicon_id[1] not in amplicons:
                amplicons[amplicon_id[1]] = {'forward':[], 'reverse':[]}
            if 'LEFT' in amplicon_id[2]:
                amplicons[amplicon_id[1]]['forward'].append(line[-1].strip().lower())
            else:
                amplicons[amplicon_id[1]]['reverse'].append(reverse_complement(line[-1].strip().lower()))
    return amplicons

def read_primerfile(filename, amplicons):
    '''
    Function that reads the primerfile of an amplicon finding run and returns the amplicons and their corresponding primersets.

    Parameters
    ----------
    filename : str
        Absolute path to the primer file.
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the ofrm of a tuple (start, end), and a dictionary containing the keys forward and reverse to store the corresponding
        forward and reverse primers.

    Returns
    -------
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the form of a tuple (start, end), and a dictionary containing the keys forward and reverse which
        contain for every amplicon the corresponding forward and reverse primers.
    primerlist : dict
        Dictionary with keys forward and reverse which are simply aggregates of the forward and reverse primers of every amplicon in $amplicons.

    '''
    primerlist = {'forward': [], 'reverse': []}
    
    with open(filename, 'r') as fn:
        lines = fn.readlines()
        line_index = 0
        while line_index < len(lines):
            cur_amplicon = lines[line_index]
            if int(cur_amplicon.split('_')[1]) > len(amplicons):
                break
            if cur_amplicon.split('_')[-1][0] == 'F': #current line corresponds to forward primer
                amplicons[int(cur_amplicon.split('_')[1]) - 1][1]['forward'].append(lines[line_index+1].strip()) #assign forward primer to corresponding amplicon
                primerlist['forward'].append(lines[line_index+1].strip()) #add primer to full set of forward primers
            else: #current line corresponds to a reverse primer
                amplicons[int(cur_amplicon.split('_')[1]) - 1][1]['reverse'].append(lines[line_index+1].strip()) #assign reverse primer to corresponding amplicon
                primerlist['reverse'].append(lines[line_index+1].strip()) #add primer to full set of reverse primers
            line_index += 2 #skip over to next amplicon
        
    return amplicons, primerlist

def locate_primers(sequence, primerlist, comparison_matrix, max_degen=10):
    '''
    Function that determines the primer binding sites of given primers in the given sequence by considering every exact binding location while
    allowing for degenerate nucleotides (i.e. if the sequence contains an N, any nucleotide will bind here).

    Parameters
    ----------
    sequence : Sequence
        Sequence object for which to determine primer binding sites.
    primerlist : dict{ 'forward':[], 'reverse':[] }
        Dictionary containing a forward and reverse key along with primers (strings) for which binding sites have to be determined.
    comparison_matrix : dict{ char : [] }
        Dictionary which dictates which nucleotides are considered equal.
    max_degen : int, optional
        Number of degenerate nucleotides allowed when checking if a primer would bind. The default is 10 (log2(4**5)).

    Returns
    -------
    hits_fwd : dict{ str : set() }
        Dictionary with the binding sites for every forward primer in $primerlist.
    hits_rev : dict{ str : set() }
        Dictionary with the binding sites for every reverse primer in $primerlist.
    '''
    sequence_reverse = reverse_complement(sequence, rev=True)
    num_degen = len(sequence)-sequence.count('a')-sequence.count('c')-sequence.count('g')-sequence.count('t')
    
    def find_hits(sequence, primers, num_degen, comparison_matrix, reverse=False, max_degen=10):
        '''
        Function that performs degenerate elastic string matching by doing k-mer based matching.

        Parameters
        ----------
        sequence : Sequence
            Sequence object for which to determine primer binding sites.
        primers : list[ str ]
            List with primers for which binding sites have to be determined.
        num_degen : int
            Number of degenerate nucleotides in the sequence. This guides the matching approach.
        comparison_matrix : dict{ char : [] }
            Dictionary which dictates which nucleotides are considered equal.
        reverse : bool, optional
            Boolean which should be set to true if the primers are based on the reverse complement of the sequence. The default is False.
        max_degen : int
            Number of degenerate nucleotides allowed in a primer binding site. The default is 10.

        Returns
        -------
        hits : dict{ str : set() }
            Dictionary with the binding sites for every primer in $primers.

        '''
        hits = {primer: set() for primer in primers}
        equivalent_characters = {'a' : ['a','r','m','w','d','h','v','n'],
                                 'c' : ['c','y','m','s','b','h','v','n'],
                                 'g' : ['g','r','k','s','b','d','v','n'],
                                 't' : ['t','y','k','w','b','d','h','n']}
        sequence_length = len(sequence)
        
        #Behaviour when there are NO degenerate nucleotides
        if num_degen == 0:
            for primer in primers:
                cont = True
                cur_occurrence = -1
                primer_length = len(primer)
                while cont:
                    try:
                        cur_occurrence = sequence.index(primer, cur_occurrence+1)
                        if not reverse:
                            hits[primer].add(cur_occurrence)
                        else:
                            hits[primer].add(sequence_length - cur_occurrence - primer_length)
                    except:
                        cont = False
        #If there are degenerate nucleotides perform a simple version of Boyer-Moore
        else:
            for primer in primers:
                cur_char_index = 0
                primer_length = len(primer)
                while cur_char_index <= sequence_length - primer_length: #iterate over sequence to find occurrence of current primer
                    stepsize = (0, False) #track where we should check for occurrence next
                    cur_degen_chars = 0 #track number of degenerate characters
                    match = True
                    for i in range(primer_length):
                        if sequence[cur_char_index + i] not in ['a','c','g','t']:
                            cur_degen_chars += 1
                        #If current character in sequence is equal to first primer character then store it as next starting point
                        if sequence[cur_char_index+i] in equivalent_characters[primer[0]] and not stepsize[1]:
                            stepsize = (max(1, i), True)
                        #If character does not match then break for-loop and start at next index
                        if not comparison_matrix[(primer[i], sequence[cur_char_index + i])][0]:
                            match = False
                            break
                    if match and cur_degen_chars <= max_degen:
                        if not reverse:
                            hits[primer].add(cur_char_index)
                        else:
                            hits[primer].add(sequence_length - cur_char_index + primer_length)
                    if stepsize[1]:
                        cur_char_index += stepsize[0]
                    else:
                        cur_char_index += i + 1
                        
        return hits
                        
    hits_fwd = find_hits(sequence, primerlist['forward'], num_degen, comparison_matrix, reverse=False)
    hits_rev = find_hits(sequence_reverse, primerlist['reverse'], num_degen, comparison_matrix, reverse=True)
    
    return hits_fwd, hits_rev    
"""
def locate_primers(sequence, primerlist, comparison_matrix, max_degen=10):
    '''
    Function that determines the primer binding sites of given primers in the given sequence by considering every exact binding location while
    allowing for degenerate nucleotides (i.e. if the sequence contains an N, any nucleotide will bind here).

    Parameters
    ----------
    sequence : Sequence
        Sequence object for which to determine primer binding sites.
    primerlist : dict{ 'forward':[], 'reverse':[] }
        Dictionary containing a forward and reverse key along with primers (strings) for which binding sites have to be determined.
    comparison_matrix : dict{ char : [] }
        Dictionary which dictates which nucleotides are considered equal.
    max_degen : int, optional
        Number of degenerate nucleotides in the sequence. This guides the k-mer based match finding method. 
        Note that this parameter is essentially redundant and can therefore be ignored. The default is 10.

    Returns
    -------
    hits_fwd : dict{ str : set() }
        Dictionary with the binding sites for every forward primer in $primerlist.
    hits_rev : dict{ str : set() }
        Dictionary with the binding sites for every reverse primer in $primerlist.
    '''
    sequence_reverse = reverse_complement(sequence, rev=True)
    max_degen = min(max_degen, len(sequence)-sequence.count('a')-sequence.count('c')-sequence.count('g')-sequence.count('t'))
    
    def find_hits(sequence, primers, max_degen, comparison_matrix, reverse=False):
        '''
        Function that performs degenerate elastic string matching by doing k-mer based matching.

        Parameters
        ----------
        sequence : Sequence
            Sequence object for which to determine primer binding sites.
        primers : list[ str ]
            List with primers for which binding sites have to be determined.
        max_degen : int
            Number of degenerate nucleotides in the sequence. This guides the matching approach.
        comparison_matrix : dict{ char : [] }
            Dictionary which dictates which nucleotides are considered equal.
        reverse : bool, optional
            Boolean which should be set to true if the primers are based on the reverse complement of the sequence. The default is False.

        Returns
        -------
        hits : dict{ str : set() }
            Dictionary with the binding sites for every primer in $primers.

        '''
        hits = {primer: set() for primer in primers}
        
        #Iterate over primers
        for primer in primers:
            min_stretch = math.ceil((len(primer) - max_degen) / (max_degen + 1)) #worst case scenario stretch of non-degenerate nucleotides to guide exact matching
            cur_kmers = [(primer[kmer:kmer+min_stretch], kmer) for kmer in range(len(primer) - min_stretch + 1)] #split primer into kmers based on worst case stretch
            for kmer, index in cur_kmers: #iterate over kmers and their starting index in the primer
                offset = -1
                cont = True
                while cont:
                    try:
                        offset = sequence.index(kmer, offset+1) #find next occurrence of kmer
                        hit = True
                        if offset >= index and len(sequence) >= offset + len(primer) - index: #check if kmer prefix and suffix can be appended in sequence
                            #Check if prefix in primer corresponds to prefix in sequence
                            for i in range(index):
                                if not comparison_matrix[(primer[i], sequence[offset-index+i])][0]: #not an exact match
                                    hit = False
                                    break
                            #Check if suffix in primer corresponds to suffix in sequence
                            if hit:
                                for i in range(len(primer) - min_stretch - index):
                                    if not comparison_matrix[(primer[index + min_stretch + i], sequence[offset + min_stretch + i])][0]: #not an exact match
                                        hit = False
                                        break
                            if hit: #prefix and suffix agree
                                if not reverse: #store occurrences
                                    hits[primer].add(offset - index)
                                else: #store occurrence translated to forward strand
                                    hits[primer].add(len(sequence) - (offset - index) - len(primer))
                                    
                    except: #no more occurrences of kmer in sequence
                        cont = False
        return hits
                        
    hits_fwd = find_hits(sequence, primerlist['forward'], max_degen, comparison_matrix, reverse=False)
    hits_rev = find_hits(sequence_reverse, primerlist['reverse'], max_degen, comparison_matrix, reverse=True)
    
    return hits_fwd, hits_rev
"""
             
def locate_amplicons(sequence, amplicons, comparison_matrix, primer_length=25, max_degen=10):
    '''
    Function that locates the realized amplicons based on the amplicons and corresponding primers in $amplicons in $sequence.
    Note that if the sequence has a stretch of $primer_length degenerate nucleotides, the results can potentially be uninterpretable if $max_degen
    is not configured appropriately.

    Parameters
    ----------
    sequence : Sequence
        Sequence object for which to determine primer binding sites.
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the form of a tuple (start, end), and a dictionary containing the keys forward and reverse which
        contain for every amplicon the corresponding forward and reverse primers in a list as values.
    comparison_matrix : dict{ char : [] }
        Dictionary which dictates which nucleotides are considered equal.
    primer_length : int, optional
        Length of the primer sequence. Theoretically this could be omitted, but since AmpliVar generates primers of fixed length
        it made sense to just include it as a parameter. The default is 25.
    max_degen : int
        Number of degenerate nucleotides in the sequence. This guides the matching approach. The default is 10 (log2(4**5)).

    Returns
    -------
    binding_sites : dict{ (int,int) : (int, int) }
        Dictionary with AmpliVar based amplicons as keys, and realized amplicons as values.

    '''
    binding_sites = {amplicon[0]: None for amplicon in amplicons}
    #Iterate over amplicons
    for amplicon in amplicons:
        fwd_hits, rev_hits = locate_primers(sequence, amplicons[amplicon], comparison_matrix)
        amplified = (0, 10**16, False)
        
        fwd_indices = set()
        rev_indices = set()
        for fwd in fwd_hits:
            fwd_indices = fwd_indices.union(fwd_hits[fwd])
        for rev in rev_hits:
            rev_indices = rev_indices.union(rev_hits[rev])
        fwd_indices = list(fwd_indices)
        rev_indices = list(rev_indices)
        #print(amplicon[0], fwd_indices, rev_indices)
        
        for fwd, rev in itertools.product(fwd_indices, rev_indices):
            #print(fwd, rev)
            if rev - fwd >= 0 and rev - fwd  < amplified[1] - amplified[0] - primer_length:
                #OLD: Amplicon including primers
                #amplified = (fwd, rev+primer_length, True)
                #NEW: Amplicon excluding primers -> reads and refs will not contain the primers
                amplified = (fwd+primer_length-1, rev, True)
        if amplified[2]:
            binding_sites[amplicon[0]] = amplified
            
    return binding_sites
"""        
DEPRECATED VERSION OF LOCATE AMPLICONS    
def locate_amplicons(sequence, amplicons, comparison_matrix, max_degen=10, primer_length=25):
    '''
    Function that locates the realized amplicons based on the amplicons and corresponding primers in $amplicons in $sequence.
    Note that if the sequence has a stretch of $primer_length degenerate nucleotides, the results can potentially be uninterpretable.

    Parameters
    ----------
    sequence : Sequence
        Sequence object for which to determine primer binding sites.
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the form of a tuple (start, end), and a dictionary containing the keys forward and reverse which
        contain for every amplicon the corresponding forward and reverse primers in a list as values.
    comparison_matrix : dict{ char : [] }
        Dictionary which dictates which nucleotides are considered equal.
    max_degen : int
        Number of degenerate nucleotides in the sequence. This guides the matching approach. The default is 10.
    primer_length : int, optional
        Length of the primer sequence. Theoretically this could be omitted, but since AmpliVar generates primers of fixed length
        it made sense to just include it as a parameter. The default is 25.

    Returns
    -------
    binding_sites : dict{ (int,int) : (int, int) }
        Dictionary with AmpliVar based amplicons as keys, and realized amplicons as values.

    '''
    binding_sites = {amplicon[0]: None for amplicon in amplicons}
    #Iterate over amplicons
    for amplicon in amplicons:
        fwd_hits, rev_hits = locate_primers(sequence, amplicon[1], comparison_matrix, max_degen=max_degen)
        amplified = (0, 10**4, False)
        
        fwd_indices = set()
        rev_indices = set()
        for fwd in fwd_hits:
            fwd_indices = fwd_indices.union(fwd_hits[fwd])
        for rev in rev_hits:
            rev_indices = rev_indices.union(rev_hits[rev])
        fwd_indices = list(fwd_indices)
        rev_indices = list(rev_indices)
        
        for fwd, rev in itertools.product(fwd_indices, rev_indices):
            if rev - fwd >= 0 and rev - fwd + primer_length < amplified[1] - amplified[0] + primer_length:
                amplified = (fwd, rev+primer_length, True)
        if amplified[2]:
            binding_sites[amplicon[0]] = amplified
            
    return binding_sites
"""


def count_primerbindings_amplicons(sequences_path, metadata_path, logfile_path, primerfile_path, max_degen=10, primer_length=25):
    amplicons = read_logfile(logfile_path) #generate amplicons from logfile
    amplicons, primerlist = read_primerfile(primerfile_path, amplicons) #refine amplicons and determine corresponding primers
    M = generate_opportunistic_matrix() #generate matrix used to determine which nucleotides are identical
    
    metadata = {} #Stores the lineage for every genome
    meta_works = True
    lineages = set()
    skip = -1
    try:
        for meta in csv.reader(open(metadata_path), delimiter='\t'):
            if skip == -1:
                for cur_meta in range(len(meta)):
                    if 'lineage' in meta[cur_meta].lower():
                        skip = cur_meta
                        break
            else:
                metadata['>' + meta[0]] = meta[skip]
                lineages.add(meta[skip])
    except Exception as e:
        print('Metadata file does not exist')
        meta_works = False
    
    primers_per_lineage = {lin: {} for lin in lineages}
    amplicons_per_lineage = {lin: {} for lin in lineages}
    total_per_lineage = {lin: 0 for lin in lineages}
    total_amplicon_bindings = {amplicon[0]: 0 for amplicon in amplicons}
    total_primer_bindings = {'forward': {}, 'reverse': {}}
    for primer in primerlist['forward']:
        total_primer_bindings['forward'][primer] = 0
    for primer in primerlist['reverse']:
        total_primer_bindings['reverse'][primer] = 0
    cur_sequence = []
    done = 0
    with open(sequences_path, 'r') as f:
        for line in f:
            if '>' in line:
                if len(cur_sequence) == 0:
                    cur_sequence = [line.strip(), '']
                else:
                    cur_fwd, cur_rev = locate_primers(cur_sequence[1], primerlist, M, max_degen=max_degen)
                    cur_amplicons = locate_amplicons(cur_sequence[1], amplicons, M, primer_length=primer_length, max_degen=max_degen)
                    if len(primers_per_lineage[metadata[cur_sequence[0]]]) == 0:
                        primers_per_lineage[metadata[cur_sequence[0]]]['forward'] = {primer: 0 for primer in cur_fwd}
                        primers_per_lineage[metadata[cur_sequence[0]]]['reverse'] = {primer: 0 for primer in cur_rev}
                    for primer in cur_fwd:
                        if len(cur_fwd[primer]) > 0:
                            primers_per_lineage[metadata[cur_sequence[0]]]['forward'][primer] += 1
                            total_primer_bindings['forward'][primer] += 1
                    for primer in cur_rev:
                        if len(cur_rev[primer]) > 0:
                            primers_per_lineage[metadata[cur_sequence[0]]]['reverse'][primer] += 1
                            total_primer_bindings['reverse'][primer] += 1
                    #Check which amplicons are amplified in current sequence
                    if len(amplicons_per_lineage[metadata[cur_sequence[0]]]) == 0:
                        amplicons_per_lineage[metadata[cur_sequence[0]]] = {amplicon: 0 for amplicon in cur_amplicons}
                    for amplicon in cur_amplicons:
                        if cur_amplicons[amplicon]:
                            amplicons_per_lineage[metadata[cur_sequence[0]]][amplicon] += 1
                            total_amplicon_bindings[amplicon] += 1
                    total_per_lineage[metadata[cur_sequence[0]]] += 1
                    cur_sequence = [line.strip(), '']
                done += 1
                print(done)
            else:
                cur_sequence[1] += line.strip().lower()
        #Process final sequence since looping over lines ends before processing final sequence     
        cur_fwd, cur_rev = locate_primers(cur_sequence[1], primerlist, M, max_degen=max_degen)
        cur_amplicons = locate_amplicons(cur_sequence[1], amplicons, M, primer_length=primer_length, max_degen=max_degen)
        if len(primers_per_lineage[metadata[cur_sequence[0]]]) == 0:
            primers_per_lineage[metadata[cur_sequence[0]]]['forward'] = {primer: 0 for primer in cur_fwd}
            primers_per_lineage[metadata[cur_sequence[0]]]['reverse'] = {primer: 0 for primer in cur_rev}
        for primer in cur_fwd:
            if len(cur_fwd[primer]) > 0:
                primers_per_lineage[metadata[cur_sequence[0]]]['forward'][primer] += 1
                total_primer_bindings['forward'][primer] += 1
        for primer in cur_rev:
            if len(cur_rev[primer]) > 0:
                primers_per_lineage[metadata[cur_sequence[0]]]['reverse'][primer] += 1
                total_primer_bindings['reverse'][primer] += 1
        #Check which amplicons are amplified in current sequence
        if len(amplicons_per_lineage[metadata[cur_sequence[0]]]) == 0:
            amplicons_per_lineage[metadata[cur_sequence[0]]] = {amplicon: 0 for amplicon in cur_amplicons}
        for amplicon in cur_amplicons:
            if cur_amplicons[amplicon]:
                amplicons_per_lineage[metadata[cur_sequence[0]]][amplicon] += 1
                total_amplicon_bindings[amplicon] += 1
        total_per_lineage[metadata[cur_sequence[0]]] += 1
        cur_sequence = [line.strip(), '']
        
    """
    #Post processing: counting primer bindings and amplifications for every lineage    
    for lineage in primers_per_lineage:
        if len(primers_per_lineage[lineage]) > 0:
            for primer in primers_per_lineage[lineage]['forward']:
                primers_per_lineage[lineage]['forward'][primer] /= total_per_lineage[lineage]
            for primer in primers_per_lineage[lineage]['reverse']:
                primers_per_lineage[lineage]['reverse'][primer] /= total_per_lineage[lineage]
        if len(amplicons_per_lineage[lineage]) > 0:
            for amplicon in amplicons_per_lineage[lineage]:
                amplicons_per_lineage[lineage][amplicon] /= total_per_lineage[lineage]
    """
    return primers_per_lineage, amplicons_per_lineage, total_per_lineage, total_primer_bindings, total_amplicon_bindings
                

def generate_simulationfile(sequences_path, metadata_path, bedfile_path, max_degen=10, max_amplicons=1000, primer_length=25):
    '''
    Function that generates a fasta file which has entries for every amplicon for the amplicons in $logfile_path, for every sequence in $sequences_path.
    The intended use for this file is to use with the ART read simulator in amplicon mode.

    Parameters
    ----------
    sequences_path : str
        Absolute path to the sequences fasta file.
    metadata_path : str
        Absolute path to the metadata file.
    logfile_path : str
        Absolute path to the log file of an AmpliVar run. For more info check read_logfile.
    primerfile_path : str
        Absolute path to the primers file of an AmpliVar run. For more info check read_primerfile
    max_degen : int, optional
        Number of degenerate nucleotides in the sequences. Note that this parameter is essentially redundant and can be ignored. The default is 10.
    max_amplicons: int, optional
        Number of amplicons to include. The default is 1000.
    primer_length : int, optional
        Length of the primer sequence. Theoretically this could be omitted, but since AmpliVar generates primers of fixed length
        it made sense to just include it as a parameter. The default is 25.

    Returns
    -------
    fasta_list : TYPE
        DESCRIPTION.

    '''
    sequences = generate_sequences(sequences_path, metadata_path) #read sequences
    amplicons = read_bedfile(bedfile_path) #generate amplicons from logfile
    M = generate_opportunistic_matrix() #generate matrix used to determine which nucleotides are identical
    
    fasta_list = []
    S = set()
    for sequence in sequences[:]:
        S.add(sequence.lineage)
        realized_amplicons = locate_amplicons(sequence.sequence_raw, amplicons, M)
        amplicon_index = 0
        for amplicon in realized_amplicons:
            amplicon_index += 1
            if realized_amplicons[amplicon]:
                fasta_list.append('>' + sequence.lineage + '_' + sequence.id + '_A' + str(amplicon_index))
                fasta_list.append(sequence.sequence_raw[realized_amplicons[amplicon][0]:realized_amplicons[amplicon][1]])
                #print(amplicon, realized_amplicons[amplicon][1] - realized_amplicons[amplicon][0])
            else:
                print('Amplicon', amplicon_index, 'not amplifiable in sequence', sequence.id)
    return fasta_list

def generate_kallistofile(sequences_path, metadata_path, bedfile_path, max_degen=10, max_amplicons=1000, primer_length=25):
    '''
    Function that generates a fasta file which has entries for every amplicon for the amplicons in $logfile_path, for every sequence in $sequences_path.
    The intended use for this file is to generate a Kallisto index using the "kallisto index" command. 

    Parameters
    ----------
    sequences_path : str
        Absolute path to the sequences fasta file.
    metadata_path : str
        Absolute path to the metadata file.
    logfile_path : str
        Absolute path to the log file of an AmpliVar run. For more info check read_logfile.
    primerfile_path : str
        Absolute path to the primers file of an AmpliVar run. For more info check read_primerfile
    max_degen : int, optional
        Number of degenerate nucleotides in the sequences. Note that this parameter is essentially redundant and can be ignored. The default is 10.
    max_amplicons: int, optional
        Number of amplicons to include. The default is 1000.
    primer_length : int, optional
        Length of the primer sequence. Theoretically this could be omitted, but since AmpliVar generates primers of fixed length
        it made sense to just include it as a parameter. The default is 25.

    Returns
    -------
    fasta_list : TYPE
        DESCRIPTION.

    '''
    sequences = generate_sequences(sequences_path, metadata_path) #read sequences
    amplicons = read_bedfile(bedfile_path) #generate amplicons from logfile
    M = generate_opportunistic_matrix() #generate matrix used to determine which nucleotides are identical
    
    fasta_list = []
    amplicon_index_dict = {} #stores which amplicons have actually been retrieved
    for sequence in sequences[:]:
        realized_amplicons = locate_amplicons(sequence.sequence_raw, amplicons, M)
        amplicon_index = 0
        cur_amplicon_indices = [] #store indices for this sequence
        s = ''
        fasta_list.append('>' + sequence.id + '\n')
        for amplicon in realized_amplicons:
            amplicon_index += 1
            if realized_amplicons[amplicon]:
                s += sequence.sequence_raw[realized_amplicons[amplicon][0]:realized_amplicons[amplicon][1]] +'A'*200
                cur_amplicon_indices.append(amplicon_index)
        fasta_list.append(s[:-200] + '\n')
        amplicon_index_dict[sequence.id] = cur_amplicon_indices
    return fasta_list, amplicon_index_dict
        
    
    

def main():
    parser = argparse.ArgumentParser(description='Generate input for ART and Kallisto')
    parser.add_argument('-s', '--sequences_path', type=str, help='Sequences file location', required=True)
    parser.add_argument('-m', '--metadata_path', type=str, help='Metadata file location', required=True)
    parser.add_argument('-b', '--bedfile_path', type=str, help='Bedfile location', required=True)
    parser.add_argument('-o', '--output_path', type=str, help='Output location (folder)', required=True)
    parser.add_argument('-n', '--num_amplicons', type=int, help='Maximum number of amplicons to consider')
    parser.add_argument('--art', action='store_true')
    parser.add_argument('--kallisto', action='store_true')
    
    args = parser.parse_args()
    if args.art:
        print('Starting with generating ART input file')
        ART_input = generate_simulationfile(args.sequences_path, args.metadata_path, args.bedfile_path)
        with open(args.output_path + '/ART_input_ARTIC.fasta', 'w') as f:
            for line in ART_input:
                f.write(line + '\n')
        print('Done with generating ART input file')
    if args.kallisto:
        print('Starting with generating Kallisto input file')
        Kallisto_input, amplicon_index_dict = generate_kallistofile(args.sequences_path, args.metadata_path, args.bedfile_path)
        with open(args.output_path + '/Kallisto_input_ARTIC.fasta', 'w') as f:
            for line in Kallisto_input:
                f.write(line)
        with open(args.output_path + '/amplicon_indices_ARTIC.csv', 'w') as f:
            for seq in amplicon_index_dict:
                cur_str = '>' + seq
                for index in amplicon_index_dict[seq]:
                    cur_str += ';' + str(index)
                f.write(cur_str + '\n')
        print('Done with generating Kallisto input file')

if __name__ == '__main__':
    main()
    '''
    #sequences, A, primerlist = main()
    F = generate_simulationfile('/Users/jaspervanbemmelen/Downloads/Netherlands_October_2022/sequences.fasta', 
                            '/Users/jaspervanbemmelen/Downloads/Netherlands_October_2022/metadata.tsv',
                            '/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/logfile_1.txt',
                            '/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/primers_1.txt')
    F2 = generate_kallistofile('/Users/jaspervanbemmelen/Downloads/Netherlands_September_2022/sequences.fasta', 
                            '/Users/jaspervanbemmelen/Downloads/Netherlands_September_2022/metadata.tsv',
                            '/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/logfile_1.txt',
                            '/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/primers_1.txt')
    sequences, primerlist, M, A, lins = main()
    sequence_bindings = [{amp[0]:0 for amp in A} for s in sequences]
    
    diff = {}
    index = 0
    fasta_list = []
    for sequence in sequences[:]:
        X = locate_amplicons(sequence, [A[0]], M)
        #print('Sequence:', index+1)
        for key in X:
            break
            print(sequence[X[key][0]:X[key][1]])
        if X[key]:
            fasta_list.append('>' + lins[index])
            fasta_list.append(sequence[X[key][0]:X[key][1]])
            subseq = sequence[X[key][0]:X[key][1]]
            if not subseq in diff:
                diff[subseq] = []
            diff[subseq].append(lins[index])
        index += 1 
    diff1 = {}
    S = []
    index = 0
    with open('/Users/jaspervanbemmelen/Downloads/Netherlands_September_2022/1673523207660.sequences.fasta', 'r') as f:
        cur = ''
        for line in f:
            if '>' in line:
                if cur != '' and cur.count('N') < 5:
                    S.append(cur.lower())
                cur = ''
            else:
                cur += line.strip()
        S.append(cur.lower())
    for s in S:
        X = locate_amplicons(s, [A[3]], M)
        for key in X:
            break
        if X[key]:
            subseq = s[X[key][0]:X[key][1]]
            if not subseq in diff:
                diff1[subseq] = []
            diff1[subseq].append(index)
        index += 1
        
    
    diff2 = {}
    S = []
    index = 0
    with open('/Users/jaspervanbemmelen/Downloads/Netherlands_October_2022/1673523305918.sequences.fasta', 'r') as f:
        cur = ''
        for line in f:
            if '>' in line:
                if cur != '' and cur.count('N') < 5:
                    print(cur)
                    S.append(cur.lower())
                cur = ''
            else:
                cur += line.strip()
        if cur.count('N') < 5:
            S.append(cur.lower())
    for s in S:
        X = locate_amplicons(s, [A[3]], M)
        for key in X:
            break
        if X[key]:
            subseq = s[X[key][0]:X[key][1]]
            if not subseq in diff:
                diff2[subseq] = []
            diff2[subseq].append(index)
        index += 1
    '''
    '''
    for amplicon in A[:]:
        index = 0
        for s in sequences[:]:
            index += 1
            fwd, rev = locate_primers(s, amplicon[1], M)
            start = 0
            start_found = False
            end = len(s)
            end_found = False
            for primer in fwd:
                if len(fwd[primer]) > 0:
                    start = max(start, max(list(fwd[primer])))
                    start_found = True
            for primer in rev:
                if len(rev[primer]) > 0:
                    end = min(end, min(list(rev[primer])))
                    end_found = True
            if start_found and end_found:
                sequence_bindings[index-1][amplicon[0]] = (start, end)
                continue
                print('Amplicon: ', amplicon)
                print(f'Realized at {start}-{end} with length {end-start}')
            else:
                print(f'For sequence {index} the amplicon {amplicon[0]} cannot be amplified')
    '''
            
    
    
    
    
    
    