from Fetch_Amplicons import generate_sequences
import argparse



def main():
    parser = argparse.ArgumentParser(description='Generates plots comparing estimation outcomes with real abundances')
    parser.add_argument('-s', '--sequences_path', type=str, help='Sequences used for simulation file location', required=True)
    parser.add_argument('-m', '--metadata_path', type=str, help='Metadata for sequences used in simulation file location', required=True)
    parser.add_argument('-a', '--abundances_path', type=str, help='TSV file containing the estimated abundances', required=True)
    parser.add_argument('-o', '--output_folder', type=str, help='Path to the folder where output will be saved', default='.')
    
    args = parser.parse_args()
    sequences = generate_sequences(args.sequences_path, args.metadata_path)
    
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
                if line[0].strip() not in estimated_abundances:
                    estimated_abundances[line[0].strip()] = float(line[-1].strip())
                    
    #Store lineages that are both in reference set and simulation set
    intersected_lineages = list(simset_lineages.intersection(refset_lineages))
    
    #Calculate difference between estimated and real abundances, and other related statistics
    MSE = 0
    MAE = 0
    for lineage in all_lineages:
        if lineage not in simset_lineages:
            errors[lineage] = estimated_abundances[lineage]
        elif lineage not in refset_lineages:
            errors[lineage] = -real_abundances[lineage]
        else:
            errors[lineage] = estimated_abundances[lineage] - real_abundances[lineage]
        super_lineage = lineage.split('.')
        if len(super_lineage) > 1:
            super_lineage = super_lineage[0] + '.' + super_lineage[1]
        else:
            super_lineage = super_lineage[0]
        super_lineages.add(super_lineage)
        if super_lineage not in super_errors:
            super_errors[super_lineage] = 0
        super_errors[super_lineage] += errors[lineage]
        MSE += (errors[lineage]**2)/len(all_lineages)
        MAE += abs(errors[lineage])/len(all_lineages)
    MSE_super = 0
    MAE_super = 0
    for lineage in super_lineages:
        MSE_super += (super_errors[lineage]**2)/len(super_lineages)
        MAE_super += abs(super_errors[lineage])/len(super_lineages)
    print('MSE:', MSE)
    print('MAE:', MAE)
    print('MSE (super)', MSE_super)
        
if __name__ == '__main__':
    main()