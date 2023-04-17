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
    
    args = parser.parse_args()
    
    max_depth=13 #this is hard-coded, change when needed
    
    for depth in range(1, max_depth+1):
        lineages = []
        errors_wgs = {}
        errors_amp = {}
        
        max_error = 0 #Used to determine xlims
        
        #Determine both intersection of lineages, and all lineages to decide which to use
        with open(args.amplicon_input + '/Seed_1/results/intersected_lineages_depth=' + str(depth) + '.csv', 'r') as f:
            for line in f:
                line = line.strip()
                if line != '':
                    lineages.append(line)
                    errors_wgs[line] = np.zeros((args.num_seeds))
                    errors_amp[line] = np.zeros((args.num_seeds))
                    
        lineages.sort(reverse=True)
                    
        for seed in range(1, args.num_seeds+1):
            #WGS
            with open(args.wgs_input + '/Seed_' + str(seed) + '/results/estimation_errors_depth=' + str(depth) + '.csv', 'r') as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        line = line.split(';')
                        if line[0] in lineages:
                            errors_wgs[line[0]][seed-1] = float(line[1])
                            max_error = max(max_error, abs(float(line[1])))
            #Amplicon
            with open(args.amplicon_input + '/Seed_' + str(seed) + '/results/estimation_errors_depth=' + str(depth) + '.csv', 'r') as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        line = line.split(';')
                        if line[0] in lineages:
                            errors_amp[line[0]][seed-1] = float(line[1])
                            max_error = max(max_error, abs(float(line[1])))
                            
        fig = plt.figure()
        ax = plt.axes()
        
        data_amp = [list(errors_amp[lineage]) for lineage in lineages]
        data_wgs = [list(errors_wgs[lineage]) for lineage in lineages]
        
        bp_left = plt.boxplot(data_wgs, positions=np.array(range(len(data_amp)))*3.0-0.4, sym='', widths=0.6, vert=False)
        bp_right = plt.boxplot(data_amp, positions=np.array(range(len(data_wgs)))*3.0+0.4, sym='', widths=0.6, vert=False)
        set_box_color(bp_left, 'red')
        set_box_color(bp_right, 'blue')
        
        plt.plot([], color='red', label='WGS')
        plt.plot([], color='blue', label='AMP')
        plt.legend()
        
        plt.yticks(range(0, len(lineages)*3, 3), lineages)
        plt.ylim(-3, len(lineages)*3)
        plt.xlim(-max_error-1, max_error+1)
        plt.tight_layout()
        
        plt.grid(color='0.8')
        plt.axvline(x=0, color='black', alpha=0.75, linewidth=0.5)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        
        plt.savefig(args.output + '/estimation_errors_plot_depth=' + str(depth) + '.png', dpi=400)

if __name__ == '__main__':
    main()