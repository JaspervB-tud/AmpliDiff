from Fetch_Amplicons import generate_sequences
import argparse



def main():
    parser = argparse.ArgumentParser(description='Generates plots comparing estimation outcomes with real abundances')
    parser.add_argument('-s', '--sequences_path', type=str, help='Sequences used for simulation file location', required=True)
    parser.add_argument('-m', '--metadata_path', type=str, help='Metadata for sequences used in simulation file location', required=True)
    parser.add_argument('-a', '--abundances_path', type=str, help='TSV file containing the estimated abundances', required=True)
    
    args = parser.parse_args()
    sequences = generate_sequences(args.sequences_path, args.metadata_path)
    
    #Store real and estimated abundances
    real_abundances = {}
    estimated_abundances = {}
    #Store lineages occurring in both simulation data (i.e. the sequences) and in the reference data (based on estimates)
    lineage_mapping = {} #just an index based map of the lineages (e.g. lineage X -> 0)
    all_lineages = set()
    simset_lineages = set([sequence.lineage for sequence in sequences])
    refset_lineages = set()
    #Read estimated abundances
    with open(args.abundances_path, 'r') as f:
        for line in f:
            line = line.split('\t')
            print(line)
    
if __name__ == '__main__':
    main()