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
    intersected_lineages = set()
    
    for lineage in estimated_abundances:
        current_lineage = '.'.join(aliases[lineage].split('.')[:depth])
        if current_lineage not in errors:
            errors[current_lineage] = 0
        if lineage not in real_abundances:
            errors[current_lineage] += estimated_abundances[lineage]
        else:
            errors[current_lineage] += estimated_abundances[lineage] - real_abundances[lineage]
            intersected_lineages.add(current_lineage)
    for lineage in real_abundances:
        current_lineage = '.'.join(aliases[lineage].split('.')[:depth])
        if current_lineage not in errors:
            errors[current_lineage] = 0
        if lineage not in estimated_abundances:
            errors[current_lineage]  += -real_abundances[lineage]
        else:
            errors[current_lineage] += estimated_abundances[lineage] - real_abundances[lineage]
            intersected_lineages.add(current_lineage)
            
    S = 0
    for lineage in errors:
        S += errors[lineage]
    print(S)
    
    return errors, intersected_lineages

"""
def determine_superlineages(lineages):
    super_lineages = {}
    for lineage in lineages:
        if len(lineage.split('.')) > 1:
            super_lineages[lineage] = lineage.split('.')[0] + '.' + lineage.split('.')[1]
        else:
            super_lineages[lineage] = lineage
    return super_lineages

def calculate_errors(lineages, super_lineages, estimated_abundances, real_abundances):
    MSE = 0
    MAE = 0
    MSE_super = 0
    MAE_super = 0
    errors = {}
    super_errors = {}
    
    for lineage in lineages:
        if lineage not in estimated_abundances:
            errors[lineage] = -real_abundances[lineage]
        elif lineage not in real_abundances:
            errors[lineage] = estimated_abundances[lineage]
        else:
            errors[lineage] = estimated_abundances[lineage] - real_abundances[lineage]
        
        super_lineage = lineage.split('.')
        if len(super_lineage) > 1:
            super_lineage = super_lineage[0] + '.' + super_lineage[1]
        else:
            super_lineage = super_lineage[0]
            
        if super_lineage in super_errors:
            super_errors[super_lineage] += errors[lineage]
        else:
            super_errors[super_lineage] = errors[lineage]
    
    for lineage in lineages:
        MSE += (errors[lineage]**2)/len(lineages)
        MAE += abs(errors[lineage])/len(lineages)
    for super_lineage in super_lineages:
        MSE_super += (super_errors[super_lineage]**2)/len(super_lineages)
        MAE_super += abs(super_errors[super_lineage])/len(super_lineages)
    
    return errors, super_errors, MSE, MSE_super, MAE, MAE_super
            
"""   

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
    TOT = 0
    for sequence in sequences:
        if sequence.lineage not in real_abundances:
            real_abundances[sequence.lineage] = 0
        real_abundances[sequence.lineage] += 100/len(sequences) #percentage of sequences per lineage
        TOT += 100/len(sequences)
        simulated_lineages.add(sequence.lineage)
        max_depth = max(max_depth, len(aliases[sequence.lineage].split('.')))
    print('total real:', TOT)
    #Determine estimated abundances
    TOT = 0
    with open(args.abundances_path, 'r') as f:
        for line in f:
            if not '#' in line:
                line = line.split('\t')
                estimated_abundances[line[0].strip()] = float(line[-1].strip())
                TOT += float(line[-1].strip())
                estimated_lineages.add(line[0].strip())
                max_depth = max(max_depth, len(aliases[line[0].strip()].split('.')))
    print('total estimated:', TOT)
    
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
    
    
    """
    #Store real and estimated abundances, and errors
    real_abundances = {} #at lineage level
    real_super_abundances = {} #at super lineage level
    estimated_abundances = {} #at lineage level
    estimated_super_abundances = {} #at super lineage level
    errors = {} #difference between estimated and actual abundance
    super_errors = {} #same as errors, but now at super lineage level (e.g. BA.1 instead of BA.1.1)
    
    #Store lineages occurring in simulation data (i.e. the sequences) or in the reference data (based on estimates)
    lineage_mapping = {} #just an index based map of the lineages (e.g. lineage X -> 0)
    all_lineages = set()
    super_lineages = set()
    simset_lineages = set()
    refset_lineages = set()
    
    #Determine actual abundances
    for sequence in sequences:
        all_lineages.add(sequence.lineage)
        simset_lineages.add(sequence.lineage)
        super_lineage = sequence.lineage.split('.')
        if len(super_lineage) > 1:
            super_lineages.add(super_lineage[0] + '.' + super_lineage[1])
        else:
            super_lineages.add(super_lineage[0])
        if sequence.lineage not in real_abundances:
            real_abundances[sequence.lineage] = 0
        real_abundances[sequence.lineage] += 100/len(sequences) #this is the percentage
        
    #Read estimated abundances
    with open(args.abundances_path, 'r') as f:
        for line in f:
            if not '#' in line:
                line = line.split('\t')
                all_lineages.add(line[0].strip())
                refset_lineages.add(line[0].strip())
                super_lineage = line[0].strip().split('.')
                if len(super_lineage) > 1:
                    super_lineages.add(super_lineage[0] + '.' + super_lineage[1])
                else:
                    super_lineages.add(super_lineage[0])
                if line[0].strip() not in estimated_abundances:
                    estimated_abundances[line[0].strip()] = float(line[-1].strip())
                    
    #Store lineages that are both in reference set and simulation set
    intersected_lineages = list(simset_lineages.intersection(refset_lineages))
    #Calculate errors and error statistics
    errors, super_errors, MSE, MSE_super, MAE, MAE_super = calculate_errors(all_lineages, super_lineages, estimated_abundances, real_abundances)

    with open(args.output_folder + '/estimation_errors.csv', 'w') as f:
        for lineage in all_lineages:
            f.write(lineage + ';' + str(errors[lineage]) + '\n')
    with open(args.output_folder + '/estimation_errors_super.csv', 'w') as f:
        for lineage in super_lineages:
            f.write(lineage + ';' + str(super_errors[lineage]) + '\n')
    with open(args.output_folder + '/intersected_lineages.txt', 'w') as f:
        for lineage in intersected_lineages:
            f.write(lineage + '\n')

    
    print('MSE (new)', MSE)
    print('MAE (new)', MAE)
    print('MSE (super) (new)', MSE_super)
    print('MAE (super) (new)', MAE_super)
    print(super_errors)
    """
if __name__ == '__main__':
    main()