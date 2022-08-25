import matplotlib.pyplot as plt
import numpy as np
from Scripts import generate_sequences

def parse_amplicons(filenames):
    amplicons = []
    for file in filenames:
        try:
            with open(file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'Amplicon' in line:
                        amplicons.append(eval(line.split(':')[-1]))
        except:
            print(file + ' does not exist')
    return amplicons

def generate_filenames(base_folder, settings, modes, parameters, parameter_names):
    res = []
    for setting in settings:
        for mode in modes:
            res.append(base_folder + setting + mode)
            
            
    progress = 0
    for parameter in parameter_names:
        cur_num_folders = len(res)
        progress = 0
        while progress < cur_num_folders:
            cur_file_name = res.pop(0)
            for variant in parameters[parameter]:
                res.append(cur_file_name + parameter + str(variant) + '_')
            progress += 1
           
    num_folders = len(res)
    progress = 0
    while progress < num_folders:
        cur_folder = res.pop(0)[:-1]
        for seed in range(1,11):
            res.append(cur_folder + '/runtimes_' + str(seed) + '.txt')
        progress += 1
    
    return res
            

if __name__ == '__main__':
    base_folder = '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/'
    settings = ['time_experiments/']
    modes = ['Soloplex/']
    parameters = {
            'ampwidth' : [200],
            'ampthresh' : [1],
            'misthresh' : [10],
            'primwidth' : [25],
            'searchwidth' : [50],
            'cov' : [1],
            'amps' : ['8_all'],
            'nseqs' : [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        }
    parameter_names = ['ampwidth', 'ampthresh', 'misthresh', 'primwidth', 'searchwidth', 'cov', 'amps', 'nseqs']
    
    filenames = generate_filenames(base_folder, settings, modes, parameters, parameter_names)
    sequences = generate_sequences('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000', 
                                   '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000')
    amplicons = parse_amplicons(filenames)
    print(amplicons)
    counts = np.zeros((sequences[0].length))
    
    for i in range(sequences[0].length):
        for a in amplicons:
            if a[0] <= i and a[1] >= i:
                counts[i] += 1
       
    fig = plt.figure(figsize=[20,10], dpi=200)
    ax = plt.gca()
    plt.title('Number of times a nucleotide is covered by an amplicon (' + str(len(filenames)) + ' runs)', size=25)
    plt.xlabel('Nucleotide index', size=18)
    plt.ylabel('Occurrence', size=18)
    plt.plot(counts, color='black', linewidth=2)
    ORF1a = [266, 13483]
    plt.axvspan(ORF1a[0], ORF1a[1], color='red', alpha=0.2)
    plt.annotate('ORF1a', ((ORF1a[0] + ORF1a[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    ORF1b = [13483, 21555]
    plt.axvspan(ORF1b[0], ORF1b[1], color='blue', alpha=0.2)
    plt.annotate('ORF1b', ((ORF1b[0] + ORF1b[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    spike_1 = [21599, 23618]
    plt.axvspan(spike_1[0], spike_1[1], color='red', alpha=0.2)
    plt.annotate('S1', ((spike_1[0]+spike_1[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    spike_2 = [23618, 25381]
    plt.axvspan(spike_2[0], spike_2[1], color='blue', alpha=0.2)
    plt.annotate('S2', ((spike_2[0]+spike_2[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    E = [26245, 26472]
    plt.axvspan(E[0], E[1], color='red', alpha=0.2)
    plt.annotate('E', ((E[0] + E[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    M = [26523, 27191]
    plt.axvspan(M[0], M[1], color='blue', alpha=0.2)
    plt.annotate('M', ((M[0] + M[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    N = [28274, 29553]
    plt.axvspan(N[0], N[1], color='red', alpha=0.2)
    plt.annotate('N', ((N[0] + N[1])/2, max(counts) + 10), color='black', alpha=0.6, size=20)
    
    y_max = max(counts) + 20
    plt.ylim([0, y_max])
    plt.savefig('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/time_experiments/Soloplex/amplicon_spread.pdf', figsize=[20,10], dpi=200, format='pdf')