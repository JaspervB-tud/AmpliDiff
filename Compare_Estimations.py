from Fetch_Amplicons import generate_sequences
import argparse

def extract_lineages(PANGO_lineages_file):
    aliases = {}
    with open(PANGO_lineages_file, 'r') as f:
        next(f)
        for line in f:
            if line[0] != '*': #if the line starts with an asterisk, the lineage has been withdrawn
                cur_line = line.split()
                if 'alias' in line or 'Alias' in line: #if the line contains alias, then the current lineage is an alias
                    aliases[cur_line[0]] = cur_line[3].replace(',', '')
                else:
                    aliases[cur_line[0]] = cur_line[0]
    return aliases

def calculate_errors(estimated_abundances, real_abundances, aliases, depth=1000):
    errors = {}
    lineages = set()
    intersected_lineages = set()
    
    for lineage in estimated_abundances:
        current_lineage = '.'.join(aliases[lineage].split('.')[:depth])
        if current_lineage not in errors:
            errors[current_lineage] = 0
        errors[current_lineage] += estimated_abundances[lineage]
        lineages.add(current_lineage)
    for lineage in real_abundances:
        current_lineage = '.'.join(aliases[lineage].split('.')[:depth])
        if current_lineage not in errors:
            errors[current_lineage] = 0
        errors[current_lineage] -= real_abundances[lineage]  
        if current_lineage in lineages:
            intersected_lineages.add(current_lineage)
    return errors, intersected_lineages

def main():
    parser = argparse.ArgumentParser(description='Generates .txt files with estimation outcome comparisons with real abundances')
    parser.add_argument('-s', '--sequences_path', type=str, help='Sequences used for simulation file location', required=True)
    parser.add_argument('-m', '--metadata_path', type=str, help='Metadata for sequences used in simulation file location', required=True)
    parser.add_argument('-a', '--abundances_path', type=str, help='TSV file containing the estimated abundances', required=True)
    parser.add_argument('--pango', type=str, help='File containing the PANGO lineage information', required=True)
    parser.add_argument('-o', '--output_folder', type=str, help='Path to the folder where output will be saved', default='.')
    
    args = parser.parse_args()
    sequences = generate_sequences(args.sequences_path, args.metadata_path) #store simulated sequences (abundance equal to occurrence)
    aliases = extract_lineages(args.pango) #store lineages and their translations
    
    real_abundances = {} #lineage level simulated abundances
    estimated_abundances = {} #lineage level extimated abundances
    
    simulated_lineages = set() #all lineages in the simulated dataset
    estimated_lineages = set() #all lineages in the reference set
    
    max_depth = 0 #maximum "length" of a lineage (e.g. B.1.1.117 has depth 4)
    
    #Determine simulated abundances
    for sequence in sequences:
        if sequence.lineage not in real_abundances:
            real_abundances[sequence.lineage] = 0
        real_abundances[sequence.lineage] += 100/len(sequences) #percentage of sequences per lineage
        simulated_lineages.add(sequence.lineage)
        max_depth = max(max_depth, len(aliases[sequence.lineage].split('.')))
    #Determine estimated abundances
    with open(args.abundances_path, 'r') as f:
        for line in f:
            if not '#' in line:
                line = line.split('\t')
                estimated_abundances[line[0].strip()] = float(line[-1].strip())
                estimated_lineages.add(line[0].strip())
                max_depth = max(max_depth, len(aliases[line[0].strip()].split('.')))
    
    #Iterate over different depths, calculate errors and store them
    for depth in range(1, max_depth+1):
        errors, cur_intersected_lineages = calculate_errors(estimated_abundances, real_abundances, aliases, depth=depth)
        MSE = 0
        MAE = 0
        print('Current depth:', depth)
        for lineage in errors:
            MSE += (errors[lineage]**2) / len(list(errors.keys()))
            MAE += abs(errors[lineage]) / len(list(errors.keys()))
        print('MSE', MSE)
        print('MAE', MAE)
        #Store estimation errors at current depth
        with open(args.output_folder + '/estimation_errors_depth=' + str(depth) + '.csv', 'w') as f:
            for lineage in errors:
                f.write(lineage + ';' + str(errors[lineage]) + '\n')
        #Store lineages in ref and sim set at current depth
        with open(args.output_folder + '/intersected_lineages_depth=' + str(depth) + '.csv', 'w') as f:
            for lineage in cur_intersected_lineages:
                f.write(lineage + '\n')

if __name__ == '__main__':
    main()