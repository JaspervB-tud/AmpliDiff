import argparse

def extract_top_n_kallisto(input_file, output_loc, n):
    content = {}
    with open(input_file, 'r') as f:
        cur_sequence = ''
        for line in f:
            if '>' in line:
                cur_sequence = line
                content[cur_sequence] = ''
            else:
                cur_amplicons = line.split('A'*200)
                for i in range(min(n, len(cur_amplicons))):
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
        cur_sequence = ''
        for line in f:
            if '>' in line:
                cur_amplicon = int(line.split('/')[-1].split('A')[-1].strip())
                if cur_amplicon <= n:
                    content.append(line.strip())
                    print(line.strip())
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
    parser.add_argument('--kallito_output', type=str, help='Output file location for kallisto')
    
    args = parser.parse_args()
    
    if args.num_amplicons > 0:
        if args.art_input != '' and args.art_output != '':
            extract_top_n_ART(args.art_input, args.art_output, args.num_amplicons)
        if args.kallisto_input != '' and args.kallisto_output != '':
            extract_top_n_kallisto(args.kallisto_input, args.kallisto_output, args.num_amplicons)

if __name__ == '__main__':
    main()