import numpy as np
from scipy import stats
import argparse

def main():
    parser = argparse.ArgumentParser(description='Generate summary statistics for all simulation results')
    parser.add_argument('-i', '--input', type=str, help='Folder containing all the results (simulation folder)', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output folder', required=True)
    
    args = parser.parse_args()
    
    amplicon_widths = ['200', '400']
    coverages = ['0.900', '0.925', '0.950', '0.975', '0.999', '1.000']
    num_amplicons = ['1', '2', '5', '10']
    depths = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']
    
    MSE = {}
    MAE = {}
    
    for coverage in coverages:
        for amplicon_width in amplicon_widths:
            for amplicons in num_amplicons:
                cur_path = args.input + '/Amplicon_' + amplicon_width + '/Coverage_' + coverage
                if coverage != '1.000':
                    cur_path += '/beta_0.05'
                cur_path += '/Simulation_' + amplicons
                for depth in depths:
                    #Calculate MAE stats
                    cur_values = []
                    with open(cur_path + '/MAE_depth=' + depth + '.tsv', 'r') as f:
                        for line in f:
                            line = line.strip()
                            if line != '':
                                cur_values.append(float(line.strip()))
                    cur_values = np.array(cur_values)
                    MAE[(coverage, amplicon_width, amplicons, depth)] = (np.mean(cur_values), np.std(cur_values, ddof=1))
                    #Calculate MSE stats
                    cur_values = []
                    with open(cur_path + '/MSE_depth=' + depth + '.tsv', 'r') as f:
                        for line in f:
                            line = line.strip()
                            if line != '':
                                cur_values.append(float(line.strip()))
                    cur_values = np.array(cur_values)
                    MSE[(coverage, amplicon_width, amplicons, depth)] = (np.mean(cur_values), np.std(cur_values, ddof=1))
                    
    header = 'Depth;'
    for coverage in coverages:
        for amplicon_width in amplicon_widths:
            for amplicon in num_amplicons:
                header += 'cov-' + coverage + ' width-' + amplicon_width + ' namps-' + amplicon + ';'
    with open(args.output + '/MSE_stats.csv', 'w') as f:
        f.write(header[:-1] + '\n')
        for depth in depths:
            cur_line = '' + depth
            for coverage in coverages:
                for amplicon_width in amplicon_widths:
                    for amplicon in num_amplicons:
                        cur_line += 'mean=' + str(MSE[(coverage, amplicon_width, amplicons, depth)]) + ', std=' + str(MSE[(coverage, amplicon_width, amplicons, depth)]) + ';'
            f.write(cur_line[:-1] + '\n')

if __name__ == '__main__':
    main()