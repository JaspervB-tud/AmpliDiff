from Fetch_Amplicons import generate_sequences
import argparse
import numpy as np

def read_errors(error_file):
    errors = {}
    with open(error_file, 'r') as f:
        for line in f:
            if ';' in line:
                cur_line = line.split(';')
                errors[cur_line[0].strip()] = float(cur_line[1].strip())
    return errors

def calculate_statistics(errors):
    MSE = np.zeros((len(errors)))
    MAE = np.zeros((len(errors)))
    
    for i, error in enumerate(errors):
        for lineage in error:
            MSE[i] += (error[lineage]**2)/len(error)
            MAE[i] += abs(error[lineage])/len(error)
            
    return MSE, MAE
        

def main():
    parser = argparse.ArgumentParser(description='Generates output based on error files in input folder')
    parser.add_argument('-i', '--input', type=str, help='Folder containing the seed folders with error files', required=True)
    parser.add_argument('-o', '--output', type=str, help='Folder where output should be stored', required=True)
    parser.add_argument('--intersect', type=str, help='File containing the lineages that appear in both ref and sim set')
    
    args = parser.parse_args()
    
    max_depth = 13 #this is hard-coded, should be adapted to particular usage
    max_seed = 20 #this is hard-coded, should be adapted to particular usage
    
    #Iterate over different depths (granularities)
    for depth in range(1,max_depth+1):
        errors = []
        #For every seed, read the errors made and calculated corresponding statistics
        for seed in range(1, max_seed+1):
            cur_error = read_errors(args.input + '/Seed_' + str(seed) + '/results/estimation_errors_depth=' + str(depth) + '.csv')
            errors.append(cur_error)
        MSE, MAE = calculate_statistics(errors)
        
        with open(args.output + '/MSE_depth=' + str(depth) + '.tsv', 'w') as f:
            for mse in MSE:
                f.write(str(mse) + '\n')
        with open(args.output + '/MAE_depth=' + str(depth) + '.tsv', 'w') as f:
            for mae in MAE:
                f.write(str(mae) + '\n')
    
if __name__ == '__main__':
    main()