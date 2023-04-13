import numpy as np
import matplotlib.pyplot as plt
import argparse

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

def main():
    parser = argparse.ArgumentParser(description='Generates plots based on error files in input folder')
    parser.add_argument('-a', '--amplicon_input', type=str, help='Folder containing the seed folders with error files for amplicon based', required=True)
    parser.add_argument('-w', '--wgs_input', type=str, help='Folder containing the seed folders with error files for wgs based', required=True)
    parser.add_argument('-o', '--output', type=str, help='Folder where output should be stored', required=True)
    parser.add_argument('-n', '--num_seeds', type=int, help='Number of seeds to include', default=20)
    parser.add_argument('--intersect', type=str, help='File containing the lineages that appear in both ref and sim set')
    
    args = parser.parse_args()
    
    lineages = []
    errors_wgs = {}
    errors_super_wgs = {}
    errors_amp = {}
    errors_super_amp = {}
    
    #Initialize error dicts and fill lineages list with lineages that are both in refset and simset
    with open(args.intersect, 'r') as f:
        for line in f:
            line = line.strip()
            if line != '':
                lineages.append(line)
                errors_wgs[line] = np.zeros((args.num_seeds))
                errors_amp[line] = np.zeros((args.num_seeds))
                
    for seed in range(1,args.num_seeds+1):
        with open(args.amplicon_input + '/Seed_' + str(seed) + '/results/estimation_errors.csv', 'r') as f:
            for line in f:
                line = line.strip()
                if line != '':
                    line = line.split(';')
                    if line[0] in lineages:
                        errors_amp[line[0]][seed-1] = float(line[1])
        with open(args.wgs_input + '/Seed_' + str(seed) + '/results/estimation_errors.csv', 'r') as f:
            for line in f:
                line = line.strip()
                if line != '':
                    line = line.split(';')
                    if line[0] in lineages:
                        errors_wgs[line[0]][seed-1] = float(line[1])
                    
    print('Lineages')
    print(lineages)
    print('*'*200)
    
    print('WGS errors')
    print(errors_wgs)
    print('*'*200)
    
    print('AMP errors')
    print(errors_amp)
    
    fig = plt.figure()
    ax = plt.axes()
    
    data_amp = [list(errors_amp[lineage]) for lineage in lineages]
    data_wgs = [list(errors_wgs[lineage]) for lineage in lineages]
    
    bp_left = plt.boxplot(data_wgs, positions=np.array(range(len(data_amp)))*2.0-0.4, sym='', widths=0.6)
    bp_right = plt.boxplot(data_amp, positions=np.array(range(len(data_wgs)))*2.0+0.4, sym='', widths=0.6)
    set_box_color(bp_left, 'red')
    set_box_color(bp_right, 'blue')
    
    plt.xticks(range(0, len(lineages)*2, 2), lineages)
    plt.xlim(-2, len(lineages)*2)
    plt.ylim(-50,50)
    plt.tight_layout()
    plt.savefig('boxcompare.png')

if __name__ == '__main__':
    main()