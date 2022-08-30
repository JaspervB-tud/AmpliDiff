from Bio import AlignIO
import csv
from Sequence import *

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