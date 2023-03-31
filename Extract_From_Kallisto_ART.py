def extract_top_n(input_file, output_loc, n):
    content = {}
    with open(input_file, 'r') as f:
        cur_sequence = ''
        for line in f:
            if '>' in line:
                cur_sequence = line
                content[cur_sequence] = ''
            else:
                cur_amplicons = line.split('A'*200)
                print(cur_amplicons)