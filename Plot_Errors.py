import numpy as np
import matplotlib.pyplot as plt
import argparse

def setBoxColors(bp):
    plt.setp(bp['boxes'][0], color='blue')
    plt.setp(bp['caps'][0], color='blue')
    plt.setp(bp['caps'][1], color='blue')
    plt.setp(bp['whiskers'][0], color='blue')
    plt.setp(bp['whiskers'][1], color='blue')
    plt.setp(bp['fliers'][0], color='blue')
    plt.setp(bp['fliers'][1], color='blue')
    plt.setp(bp['medians'][0], color='blue')

    plt.setp(bp['boxes'][1], color='red')
    plt.setp(bp['caps'][2], color='red')
    plt.setp(bp['caps'][3], color='red')
    plt.setp(bp['whiskers'][2], color='red')
    plt.setp(bp['whiskers'][3], color='red')
    plt.setp(bp['fliers'][2], color='red')
    plt.setp(bp['fliers'][3], color='red')
    plt.setp(bp['medians'][1], color='red')

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
    plt.hold(True)
    for i in range(len(lineages)):
        bp = plt.boxplot([errors_wgs[lineages[i]], errors_amp[lineages[i]]], positions=[1+3*i, 2+3*i], widths=0.5)
        plt.setBoxColors(bp)
    plt.xlim(0, 3*len(lineages))
    plt.ylim(0, 100)
    ax.set_xticklabels(lineages)
    ax.set_xticks([1.5*(i+1) for i in range(len(lineages))])
    plt.savefig('boxcompare.png')

if __name__ == '__main__':
    main()