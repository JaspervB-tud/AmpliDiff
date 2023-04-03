# Actual imports
import argparse
from Class_Methods import generate_sequences
from Sequence import *
import numpy as np
import Bio
from Bio import pairwise2
from pathlib import Path
import matplotlib.pyplot as plt

################################# ACTUAL CODE #######################################################
def parse_logfile(location, num_positions):
    positions_covered = np.zeros((num_positions), dtype=np.int16)
    with open(location, 'r') as f:
        for line in f:
            if 'succesfully' in line:
                amplicon = line.split(')')[0].split('(')[1].split(',')
                amplicon = [int(amplicon[0].strip()), int(amplicon[1].strip())]
                positions_covered[amplicon[0]:amplicon[1]] += 1
    return positions_covered

#################################### HELPER FUNCTIONS ###################################################
def determine_annotations(sequences, ref_genome):
    alignments = pairwise2.align.globalxx(sequences[0].sequence_raw.upper(), str(ref_genome.seq))
    ref_sequence = Sequence(alignments[0].seqB.lower(), 'NC_045512.2')
    ref_sequence.align_to_trim()
    sequences[0].align_to_trim()
    
    #Annotations are based on annotations for reference genome "NC_045512.2"
    annotations = {}
    annotations2 = {}
    ORF1a = [266, 13483]
    ORF1b = [13483, 21555]
    spike_1 = [21599, 23618]
    spike_2 = [23618, 25381]
    E = [26245, 26472]
    M = [26523, 27191]
    ORF6 = [27202, 27384]
    ORF7a = [27439, 27756]
    ORF7b = [27756, 27884]
    N = [28274, 29553]
    
    #ORF1a
    annotations['ORF1a'] = [np.where(ref_sequence.aligned_to_trim == ORF1a[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF1a[1])[0][0]]
    annotations2['ORF1a'] = [np.where(sequences[0].aligned_to_trim == annotations['ORF1a'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['ORF1a'][1])[0][0]]
    #ORF1b
    annotations['ORF1b'] = [np.where(ref_sequence.aligned_to_trim == ORF1b[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF1b[1])[0][0]]
    annotations2['ORF1b'] = [np.where(sequences[0].aligned_to_trim == annotations['ORF1b'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['ORF1b'][1])[0][0]]
    #S1
    annotations['S1'] = [np.where(ref_sequence.aligned_to_trim == spike_1[0])[0][0], np.where(ref_sequence.aligned_to_trim == spike_1[1])[0][0]]
    annotations2['S1'] = [np.where(sequences[0].aligned_to_trim == annotations['S1'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['S1'][1])[0][0]]
    #S2
    annotations['S2'] = [np.where(ref_sequence.aligned_to_trim == spike_2[0])[0][0], np.where(ref_sequence.aligned_to_trim == spike_2[1])[0][0]]
    annotations2['S2'] = [np.where(sequences[0].aligned_to_trim == annotations['S2'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['S2'][1])[0][0]]
    #E
    annotations['E'] = [np.where(ref_sequence.aligned_to_trim == E[0])[0][0], np.where(ref_sequence.aligned_to_trim == E[1])[0][0]]
    annotations2['E'] = [np.where(sequences[0].aligned_to_trim == annotations['E'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['E'][1])[0][0]]
    #M
    annotations['M'] = [np.where(ref_sequence.aligned_to_trim == M[0])[0][0], np.where(ref_sequence.aligned_to_trim == M[1])[0][0]]
    annotations2['M'] = [np.where(sequences[0].aligned_to_trim == annotations['M'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['M'][1])[0][0]]
    #ORF6
    annotations['ORF6'] = [np.where(ref_sequence.aligned_to_trim == ORF6[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF6[1])[0][0]]
    annotations2['ORF6'] = [np.where(sequences[0].aligned_to_trim == annotations['ORF6'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['ORF6'][1])[0][0]]
    #ORF7a
    annotations['ORF7a'] = [np.where(ref_sequence.aligned_to_trim == ORF7a[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF7a[1])[0][0]]
    annotations2['ORF7a'] = [np.where(sequences[0].aligned_to_trim == annotations['ORF7a'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['ORF7a'][1])[0][0]]
    #ORF7b
    annotations['ORF7b'] = [np.where(ref_sequence.aligned_to_trim == ORF7b[0])[0][0], np.where(ref_sequence.aligned_to_trim == ORF7b[1])[0][0]]
    annotations2['ORF7b'] = [np.where(sequences[0].aligned_to_trim == annotations['ORF7b'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['ORF7b'][1])[0][0]]
    #ORF6 + 7
    annotations2['ORF6-7'] = [annotations2['ORF6'][0], annotations2['ORF7b'][1]]
    #N
    annotations['N'] = [np.where(ref_sequence.aligned_to_trim == N[0])[0][0], np.where(ref_sequence.aligned_to_trim == N[1])[0][0]]
    annotations2['N'] = [np.where(sequences[0].aligned_to_trim == annotations['N'][0])[0][0], np.where(sequences[0].aligned_to_trim == annotations['N'][1])[0][0]]
    
    return annotations2
#####################################################################################################

def main():
    parser = argparse.ArgumentParser(description='Plot amplicon coverage for given amplicon width coverage and beta value')
    parser.add_argument('-w', '--ampwidth', type=str, required=True)
    parser.add_argument('-c', '--coverage', type=str, required=True)
    parser.add_argument('-b', '--beta', type=str, required=True)
    parser.add_argument('-o', '--output_folder', type=str, required=True)
    
    args = parser.parse_args()
    
    # Get sequences in order to align to reference genome to allow region annotation for plots
    print('Reading sequences')
    sequences_folder = '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasmijn/ref_sets_gisaid_2022_08_18/global_all_time_N0_L29000'
    sequences = generate_sequences(sequences_folder, sequences_folder)
    num_positions = sequences[0].length
    print('Done reading sequences')
    
    # Get reference genome and determine annotations by aligning
    print('Reading reference genome and determining annotations')
    ref_genome = Bio.SeqIO.read(open('/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/amplivar/NC_045512.2.fasta'), format='fasta')
    annotations = determine_annotations(sequences, ref_genome)
    regions = ['ORF1a', 'ORF1b', 'S1', 'S2', 'E', 'M', 'ORF6-7', 'N']
    colors = ['red', 'blue']
    print('Done reading reference genome and determining annotations')
    
    # Folder containing subdirectories with logfiles
    base_folder = '/tudelft.net/staff-umbrella/SARSCoV2Wastewater/jasper/source_code/final_scripts/fast_output/Global/downsampling/Soloplex/'
    positions_covered_all = np.zeros((num_positions), dtype=np.int16) #store how often each position is covered in every run
    total_runs = 0 #store the number of succesful runs
    # Iterate over number of sequences considered in downsampling runs
    for num_seqs in [1500, 2000, 2500]:
        print('Working on downsampling with', num_seqs, 'sequences')
        positions_covered_cur = np.zeros((num_positions), dtype=np.int16)
        actual_runs = 0
        for seed in range(1,16):
            if args.coverage == '1.000':
                logfile = (base_folder + 'coverage-' + args.coverage + '/amplicon_width-' + args.ampwidth +
                                                       '/primer_width-25/amplicon_threshold-1/misthresh' + args.ampwidth[:2] +
                                                       '_searchwidth50_amps10_all_nseqs' + str(num_seqs) + '/logfile_' + str(seed) + '.txt')
            else:
                logfile = (base_folder + 'coverage-' + args.coverage + '/beta-' + args.beta + '/amplicon_width-' + args.ampwidth +
                                                       '/primer_width-25/amplicon_threshold-1/misthresh' + args.ampwidth[:2] +
                                                       '_searchwidth50_amps10_all_nseqs' + str(num_seqs) + '/logfile_' + str(seed) + '.txt')
            try:
                positions_covered_cur += parse_logfile(logfile, num_positions)
                actual_runs += 1
            except Exception as e:
                print(e)
                continue
        positions_covered_all += positions_covered_cur #add covered positions for current number of sequences to aggregated total
        total_runs += actual_runs
        
        if actual_runs > 0:
            color_index = 0
            fig = plt.figure(figsize=[20,10], dpi=200)
            ax = plt.gca()
            plt.title('Coverage for amplicons of width 400 (' + str(actual_runs) + ' runs) while subsampling ' + str(num_seqs) + ' sequences', size=25)
            plt.xlabel('Nucleotide index', size=20)
            plt.ylabel('Relative coverage', size=20)
            print(positions_covered_cur)
            print('Plotting')
            plt.plot(positions_covered_cur/actual_runs, color='black', linewidth=3)
            for region in regions:
                plt.axvspan(annotations[region][0], annotations[region][1], color=colors[color_index % 2], alpha=0.2)
                plt.annotate(region, ((annotations[region][0] + annotations[region][1])/2, 0.9), color='black', alpha=0.6, size=20, ha='center', rotation=90)
                color_index += 1
            plt.ylim([0, 1.2])
            if args.coverage == '1.000':
                output_loc = (args.output_folder + '/coverage-' + args.coverage + '_ampliconwidth-' + args.ampwidth + 
                              '_primerwidth-25_ampliconthreshold-1_misthresh-' + args.ampwidth[:2] + '_searchwidth-50_amps-10_nseqs-' + 
                              str(num_seqs) + '.pdf') 
            else:
                output_loc = (args.output_folder + '/coverage-' + args.coverage + '/beta-' + args.beta + '_ampliconwidth-' + args.ampwidth + 
                              '_primerwidth-25_ampliconthreshold-1_misthresh-' + args.ampwidth[:2] + '_searchwidth-50_amps-10_nseqs-' + 
                              str(num_seqs) + '.pdf') 
            plt.savefig(output_loc, figsize=[20,10], dpi=200, format='pdf')
            del fig, ax
            
    
if __name__ == '__main__':
    main()
    

    