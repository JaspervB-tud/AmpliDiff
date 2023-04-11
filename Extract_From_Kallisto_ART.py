import argparse

def extract_top_n_kallisto(input_file, amplicon_index_file, output_loc, n):
    amplicon_indices = {}
    with open(amplicon_index_file, 'r') as f:
        for line in f:
            if '>' in line:
                cur_line = line.split(';')
                amplicon_indices[cur_line[0]] = [int(cur_line[i].strip()) for i in range(1,len(cur_line))]
    
    content = {}
    with open(input_file, 'r') as f:
        cur_sequence = ''
        for line in f:
            #Check if line is a header
            if '>' in line:
                cur_sequence = line
                content[cur_sequence] = ''
            #If line is not a header, add amplicons
            else:
                cur_amplicons = line.split('A'*200)
                for i in range(min(n, len(cur_amplicons))): #iterate over amplicons
                    if amplicon_indices[cur_sequence.strip()] <= n: #check if current amplicon should be included -> this may imply that some sequence get no amplicons
                        content[cur_sequence] += cur_amplicons[i].strip() + 'A'*200
                content[cur_sequence] = content[cur_sequence][:-200]
    s = ''
    for sequence in content:
        s += sequence
        s += content[sequence] + '\n'
    with open(output_loc, 'w') as f:
        f.write(s.rstrip())
    return content

def extract_top_n_ART(input_file, output_loc, n):
    content = []
    with open(input_file, 'r') as f:
        for line in f:
            if '>' in line:
                cur_amplicon = int(line.split('/')[-1].split('A')[-1].strip())
                if cur_amplicon <= n:
                    content.append(line.strip())
            else:
                if cur_amplicon <= n:
                    content.append(line.strip())
    s = ''
    for entry in content:
        s += entry + '\n'
    with open(output_loc, 'w') as f:
        f.write(s.rstrip())
    return content
                
def main():
    parser = argparse.ArgumentParser(description='Use ART and Kallisto input files to generate input files using only the first n amplicons')
    parser.add_argument('-n', '--num_amplicons', type=int, help='Consider the first n amplicons', default=0)
    parser.add_argument('--art_input', type=str, help='ART_input file location', default='')
    parser.add_argument('--art_output', type=str, help='Output file location for ART', default='')
    parser.add_argument('--kallisto_input', type=str, help='Kallisto_input file location')
    parser.add_argument('--amplicon_index', type=str, help='Amplicon indices that could be amplified for every reference genome (Kallisto only)', default='')
    parser.add_argument('--kallisto_output', type=str, help='Output file location for kallisto')
    
    args = parser.parse_args()
    
    if args.num_amplicons > 0:
        if args.art_input != '' and args.art_output != '':
            print('Extracting ART amplicons')
            extract_top_n_ART(args.art_input, args.art_output, args.num_amplicons)
        if args.kallisto_input != '' and args.kallisto_output != '':
            print('Extracting Kallisto amplicons')
            extract_top_n_kallisto(args.kallisto_input, args.amplicon_index, args.kallisto_output, args.num_amplicons)

if __name__ == '__main__':
    main()