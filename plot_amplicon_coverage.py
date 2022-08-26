import matplotlib.pyplot as plt
import numpy as np
from Scripts import generate_sequences
from Classes import Sequence
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import Bio

def parse_amplicons(filenames, sequence_length):
    counts_200 = np.zeros((sequence_length))
    counts_400 = np.zeros((sequence_length))
    counts = np.zeros((sequence_length))
    n_files = [0,0]
    for file in filenames:
        try:
            with open(file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'Amplicon' in line:
                        cur_amplicon = eval(line.split(':')[-1])
                        if 'ampwidth200' in file:
                            for i in range(cur_amplicon[0], cur_amplicon[1]):
                                counts_200[i] += 1
                                counts[i] += 1
                        elif 'ampwidth400' in file:
                            for i in range(cur_amplicon[0], cur_amplicon[1]):
                                counts_400[i] += 1
                                counts[i] += 1
            if 'ampwidth200' in file:
                n_files[0] += 1
            if 'ampwidth400' in file:
                n_files[1] += 1
        except:
            print(file + ' does not exist')
    print('Done determining occurrences')
    return counts_200, counts_400, counts, n_files

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

def determine_annotations(sequences, ref_genome):
    alignments = pairwise2.align.globalxx(sequences[0].sequence_raw.upper(), str(ref_genome.seq))
    ref_sequence = Sequence(alignments[0].seqB.lower(), 'NC_045512.2')
    ref_sequence.align_to_trim()
    
    annotations = {}
    ORF1a = [266, 13483]
    ORF1b = [13483, 21555]
    spike_1 = [21599, 23618]
    spike_2 = [23618, 25381]
    E = [26245, 26472]
    M = [26523, 27191]
    N = [28274, 29553]
    annotations['ORF1a'] = [np.where(ref_sequence.aligned_to_trim == ORF1a[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF1a[1])[0][0]]
    annotations['ORF1b'] = [np.where(ref_sequence.aligned_to_trim == ORF1b[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF1b[1])[0][0]]
    annotations['S1'] = [np.where(ref_sequence.aligned_to_trim == spike_1[0])[0][0], np.where(ref_sequence.aligned_to_trim == spike_1[1])[0][0]]
    annotations['S2'] = [np.where(ref_sequence.aligned_to_trim == spike_2[0])[0][0], np.where(ref_sequence.aligned_to_trim == spike_2[1])[0][0]]
    annotations['E'] = [np.where(ref_sequence.aligned_to_trim == E[0])[0][0], np.where(ref_sequence.aligned_to_trim == E[1])[0][0]]
    annotations['M'] = [np.where(ref_sequence.aligned_to_trim == M[0])[0][0], np.where(ref_sequence.aligned_to_trim == M[1])[0][0]]
    annotations['N'] = [np.where(ref_sequence.aligned_to_trim == N[0])[0][0], np.where(ref_sequence.aligned_to_trim == N[1])[0][0]]
    print(annotations)
    return annotations
            

if __name__ == '__main__':
    base_folder = '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/'
    settings = ['time_experiments/']
    modes = ['Soloplex/']
    parameters = {
            'ampwidth' : [200, 400],
            'ampthresh' : [1],
            'misthresh' : [10],
            'primwidth' : [25],
            'searchwidth' : [50],
            'cov' : [1],
            'amps' : ['8_all'],
            'nseqs' : [500, 600, 700, 800, 900, 1000]
        }
    parameter_names = ['ampwidth', 'ampthresh', 'misthresh', 'primwidth', 'searchwidth', 'cov', 'amps', 'nseqs']
    
    filenames = generate_filenames(base_folder, settings, modes, parameters, parameter_names)
    sequences = generate_sequences('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000', 
                                   '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000')
    ref_genome = Bio.SeqIO.read(open('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/amplivar/NC_045512.2.fasta'), format='fasta')
    counts_200, counts_400, counts, n_files = parse_amplicons(filenames, sequences[0].length)
                    
    annotations = determine_annotations(sequences, ref_genome)
    regions = ['ORF1a', 'ORF1b', 'S1', 'S2', 'E', 'M', 'N']
    colors = ['red','blue']
    color_index = 0
    
    #Aggregated counts 
    fig = plt.figure(figsize=[20,10], dpi=200)
    ax = plt.gca()
    plt.title('Number of times a nucleotide is covered by an amplicon (' + str(n_files[0] + n_files[1]) + ' runs)', size=25)
    plt.xlabel('Nucleotide index', size=18)
    plt.ylabel('Occurrence', size=18)
    plt.plot(counts/(n_files[0] + n_files[1]), color='black', linewidth=2)
    for region in regions:
        plt.axvspan(annotations[region][0], annotations[region][1], color=colors[color_index % 2], alpha=0.2)
        plt.annotate(region, ((annotations[region][0] + annotations[region][1])/2, 1.1, color='black', alpha=0.6, size=20, ha='center')
        color_index += 1
    plt.ylim([0, 1.2])
    plt.savefig('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/time_experiments/Soloplex/amplicon_spread_all.pdf', figsize=[20,10], dpi=200, format='pdf')
    del fig, ax
    
    
    #200-width counts
    fig = plt.figure(figsize=[20,10], dpi=200)
    ax = plt.gca()
    plt.title('Number of times a nucleotide is covered by an amplicon (' + str(n_files[0]) + ' runs) of width 200', size=25)
    plt.xlabel('Nucleotide index', size=18)
    plt.ylabel('Occurrence', size=18)
    plt.plot(counts_200/n_files[0], color='black', linewidth=2)
    for region in regions:
        plt.axvspan(annotations[region][0], annotations[region][1], color=colors[color_index % 2], alpha=0.2)
        plt.annotate(region, ((annotations[region][0] + annotations[region][1])/2, 1.1, color='black', alpha=0.6, size=20, ha='center')
        color_index += 1
    plt.ylim([0, 1.2)
    plt.savefig('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/time_experiments/Soloplex/amplicon_spread_200.pdf', figsize=[20,10], dpi=200, format='pdf')
    del fig, ax
    
    #400-width counts
    fig = plt.figure(figsize=[20,10], dpi=200)
    ax = plt.gca()
    plt.title('Number of times a nucleotide is covered by an amplicon (' + str(n_files[1]) + ' runs) of width 400', size=25)
    plt.xlabel('Nucleotide index', size=18)
    plt.ylabel('Occurrence', size=18)
    plt.plot(counts_400/n_files[1], color='black', linewidth=2)
    for region in regions:
        plt.axvspan(annotations[region][0], annotations[region][1], color=colors[color_index % 2], alpha=0.2)
        plt.annotate(region, ((annotations[region][0] + annotations[region][1])/2, 1.1), color='black', alpha=0.6, size=20, ha='center')
        color_index += 1
    plt.ylim([0, 1.2)
    plt.savefig('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/time_experiments/Soloplex/amplicon_spread_400.pdf', figsize=[20,10], dpi=200, format='pdf')
    
    