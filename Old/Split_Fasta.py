import argparse

def split_sequences(sequences_path, sequences_per_file, output_path):
    num_processed = 0
    batch = 1
    cur_sequences = ''
    with open(sequences_path, 'r') as f:
        for line in f:
            if '>' in line:
                if cur_sequences == '': #First sequence in the file
                    cur_sequences += line.split('|')[0] + '\n' #Headers are formatted as ID|DATE thus we need to remove the date to retain the actual ID
                else:
                    num_processed += 1
                    if num_processed % sequences_per_file == 0: #write to new file
                        with open(output_path + '/sequences_batch_' + str(batch) + '.fasta', 'w') as o:
                            o.write(cur_sequences)
                        batch += 1
                        cur_sequences = line.split('|')[0] + '\n'
                    else:
                        cur_sequences += line.split('|')[0] + '\n'
            else:
                cur_sequences += line
        if len(cur_sequences) > 0:
            with open(output_path + '/sequences_batch_' + str(batch) + '.fasta', 'w') as o:
                o.write(cur_sequences)
                
def main():
    parser = argparse.ArgumentParser(description='Split genomes in fasta file into different files for use in multiprocessing')
    parser.add_argument('-s', '--sequences_path', type=str, help='Sequences file location', required=True)
    parser.add_argument('-b', '--batch_size', type=int, help='Number of sequences per output file', required=True)
    parser.add_argument('-o', '--output_path', type=str, help='Output folder', required=True)
    args = parser.parse_args()
    
    split_sequences(args.sequences_path, args.batch_size, args.output_path)
    
if __name__ == '__main__':
    main()