import argparse
from Fetch_Amplicons import generate_sequences

def convert_metadata(sequences_path, metadata_path, output_path):
    '''
    Function that reformats the metadata obtained when downloading subsets of the GISAID database to the format obtained 
    when downloading the entire database.

    Parameters
    ----------
    sequences_path : str
        Absolute path to the sequences described in the metadata file
    metadata_path : str
        Absolute path to the "wrong" metadata file, assumption is that metadata has the same order as the sequence file.
    output_path : str
        Absolute path to where the output should be stored (including filename and extension).

    Returns
    -------
    None.

    '''
    correct_headers = ['Virus name', 'Type', 'Accession ID', 'Collection date', 'Location', 'Additional location information', 'Sequence length',
                   'Host', 'Patient age', 'Gender', 'Clade', 'Pango lineage', 'Pangolin version', 'Variant', 'AA Substitutions',
                   'Submission date', 'Is reference?', 'Is complete?', 'Is high coverage?', 'Is low coverage?', 'N-Content', 'GC-content',
                   'date']
    reformatted_meta = ['\t'.join(correct_headers)]
    sequences = generate_sequences(sequences_path, metadata_path, max_n=10**6)
    with open(metadata_path, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            cur_line = line.split('\t')
            cur_sequence = cur_line[0] + '\t' + cur_line[1] + '\t' + cur_line[2] + '\t' + cur_line[26] + '\t'
            cur_location = ''
            for loc_line_index in range(5,9):
                if cur_line[loc_line_index].strip() != '':
                    if cur_location == '':
                        cur_location += cur_line[loc_line_index]
                    else:
                        cur_location += ' / ' + cur_line[loc_line_index]
            #0.0, 0.0 -> N-content, GC-content (relative)
            cur_sequence += cur_location + '\t' + '\t' + cur_line[13] + '\t' + cur_line[14] + '\t' + cur_line[15] + '\t' + \
                                cur_line[16] + '\t' + cur_line[19] + '\t' + cur_line[18] + '\t' + 'unknown' + '\t' + 'insert variant here' + '\t' + \
                                'insert AA substitutions here' + '\t' + cur_line[26] + '\t' + '\t' + 'TRUE' + '\t' + 'TRUE' + '\t' + '\t' + \
                                str(sequences[sequences.index(cur_line[0])].sequence_raw.lower().count('n')/len(sequences[sequences.index(cur_line[0])].sequence_raw)) + '\t' + '0.0' + '\t' + cur_line[4]
            reformatted_meta.append(cur_sequence.split('\t'))
            #print(cur_sequence.split('\t'))
    with open(output_path, 'w') as f:
        f.write('\t'.join(correct_headers) + '\n')
        for line in reformatted_meta[1:]:
            f.write('\t'.join(line) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Convert the metadata obtained from GISAID subset to proper format required for VLQ')
    parser.add_argument('-i', '--input_path', type=str, help='Metadata file location', required=True)
    parser.add_argument('-o', '--output_path', type=str, help='Output location', required=True)
    parser.add_argument('-s', '--sequences_path', type=str, help='Sequences file location', required=True)
    args = parser.parse_args()
    convert_metadata(args.sequences_path, args.input_path, args.output_path)
    '''
    metadata_path = '/Users/jaspervanbemmelen/Documents/Wastewater/Data/full_data_set/metadata.tsv'
    with open(metadata_path, 'r') as f:
        header = f.readline().split('\t')
        print(header)
    print()
    with open('/Users/jaspervanbemmelen/Downloads/Netherlands_September_2022/metadata.tsv', 'r') as f:
        header = f.readline().split('\t')
        for i in range(len(header)):
            print(i, ':', header[i])
        print('-'*100)
        print(f.readline())
    print('-'*100)
    convert_metadata('/Users/jaspervanbemmelen/Downloads/Netherlands_September_2022/metadata.tsv', '/Users/jaspervanbemmelen/Documents/test.tsv')
    '''
if __name__ == '__main__':
    main()