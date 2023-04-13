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
    
    means = {}
    stds = {}
    
    for error in errors:
        for lineage in error:
            if len(lineage) > 0:
                if lineage not in means:
                    means[lineage] = [error[lineage]]
                else:
                    means[lineage].append(error[lineage])
    i = 0
    for error in errors:
        for lineage in error:
            MSE[i] += (error[lineage]**2)/len(error)
            MAE[i] += abs(error[lineage])/len(error)
            cur_errors = np.array(means[lineage])
            try:
                stds[lineage] = np.std(cur_errors, ddof=1)
            except:
                stds[lineage] = 0
            means[lineage] = np.mean(cur_errors)
        i += 1
    return means, stds, MSE, MAE
        

def main():
    parser = argparse.ArgumentParser(description='Generates plots based on error files in input folder')
    parser.add_argument('-i', '--input', type=str, help='Folder containing the seed folders with error files', required=True)
    parser.add_argument('-o', '--output', type=str, help='Folder where output should be stored', required=True)
    parser.add_argument('--intersect', type=str, help='File containing the lineages that appear in both ref and sim set')
    
    args = parser.parse_args()
    
    errors = []
    super_errors = []
    for seed in range(1,21):
        cur_error = read_errors(args.input + '/Seed_' + str(seed) + '/results/estimation_errors.csv')
        cur_super_error = read_errors(args.input + '/Seed_' + str(seed) + '/results/estimation_errors_super.csv')
        errors.append(cur_error)
        super_errors.append(cur_super_error)
    mu, sigma, MSE, MAE = calculate_statistics(errors)
    mu, sigma, MSE_super, MAE_super = calculate_statistics(super_errors)
    
    with open(args.output + '/MSE.tsv', 'w') as f:
        for mse in MSE:
            f.write(str(mse) + '\n')
    with open(args.output + '/MAE.tsv', 'w') as f:
        for mae in MAE:
            f.write(str(mae) + '\n')
    with open(args.output + '/MSE_super.tsv', 'w') as f:
        for mse_super in MSE_super:
            f.write(str(mse_super) + '\n')
    with open(args.output + '/MAE_super.tsv', 'w') as f:
        for mae_super in MAE_super:
            f.write(str(mae_super) + '\n')
        
    
if __name__ == '__main__':
    main()