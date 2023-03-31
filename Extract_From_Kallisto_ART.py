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
                
content3 = extract_top_n_kallisto('/Users/jaspervanbemmelen/Documents/Wastewater/Kallisto_input.fasta', '/Users/jaspervanbemmelen/Documents/Wastewater/Kallisto_test.fasta', 5)
content4 = extract_top_n_ART('/Users/jaspervanbemmelen/Documents/Wastewater/ART_input.fasta', '/Users/jaspervanbemmelen/Documents/Wastewater/ART_test.fasta', 5)
