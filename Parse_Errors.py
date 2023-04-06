from Fetch_Amplicons import generate_sequences
import argparse

def read_errors(error_file):
    errors = {}
    with open(error_file, 'r') as f:
        for line in f:
            if ';' in line:
                cur_line = line.split(';')
                errors[cur_line[0].strip()] = cur_line[1].strip()
    return errors

def main():
    parser = argparse.ArgumentParser(description='Generates plots based on error files in input folder')
    parser.add_argument('-i', '--input', type=str, help='Folder containing the seed folders with error files')
    parser.add_argument('-o', '--output', type=str, help='Folder where output should be stored')
    
    args = parser.parse_args()
    
    errors = []
    super_errors = []
    for seed in range(1,21):
        cur_error = read_errors(args.input + '/Seed_' + str(seed) + '/results/estimation_errors.csv')
        cur_super_error = read_errors(args.input + '/Seed_' + str(seed) + '/results/estimation_errors_super.csv')
        errors.append(cur_error)
        super_errors.append(cur_super_error)
    print('Errors')
    print(errors)
    print('Super errors')
    print(super_errors)
    
if __name__ == '__main__':
    main()