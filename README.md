# AmpliDiff

This repository contains the code for [AmpliDiff](https://www.biorxiv.org/content/10.1101/2023.07.22.550164v1), a Python tool that finds amplicons and corresponding primers in viral genomes in order to differentiate between different lineages/strains/species.

### Requirements
All dependencies can be found in the dependencies_AmpliDiff.txt and environment_AmpliDiff.yml files. Note that AmpliDiff has not been tested on Windows systems and it is therefore discouraged to run AmpliDiff on Windows.

### Installation
AmpliDiff can be installed by cloning this repo, installing the dependencies (through conda using the environment_AmpliDiff.yml file for example) and building the amplicon_generation.pyx Cython file using the following commands:
```
git clone https://github.com/JaspervB-tud/AmpliDiff.git
cd AmpliDiff/AmpliDiff
python setup.py build_ext --inplace
```
As of now, AmpliDiff uses [Gurobi](https://www.gurobi.com) to solve the primer feasibility and minimization problems.

### Usage
AmpliDiff requires the user to provide an Multiple Sequence Alignment in the form of a fasta file, and a metadata file in TSV format which has a header line containing a header with "lineage" that contains the lineages/strains/species of an entry. The first column should always be a sequence identifier which is used to assign classes to genomes in the MSA fasta file. Note that sequence IDs in the metadata should not contain "|" characters and should correspond exactly to the IDs in the fasta file.

AmpliDiff then has the following list of optional parameters:
```
-o                     : Path to folder where output will be stored. Default is current folder.
--primer_thresholds  : Path to the primer thresholds file. Default is ./primer_thresholds.csv.
##Amplicon parameters
-aw                    : Amplicon width. Default is 200.
-mm                    : Number of allowed mismatches during amplicon differentiation. Default is 0.
-mt                    : Number of allowed misalignment characters in an amplicon. Default is 20.
##Primer parameters
-pw                    : Primer size. Default is 25.
-sw                    : Search window flanking an amplicon. Default is 50.
-cov                   : Minimal required amplifiability. Default is 1 (100%).
-b                     : Trade-off parameter between primer pairs and differentiability. Default is 0.05.
--max_primer_degeneracy: Maximimum allowed degeneracy for disambiguating primer candidates. Default is 1024.
--gc_lb                : Minimum required GC-content in primers. Default is 0.4 (40%).
--gc_ub                : Maximum allowed GC-content in primers. Default is 0.6 (60%).
--melting_lb           : Minimum required primer melting temperature in degrees Celsius. Default is 55.
--melting_ub           : Maximum allowed primer melting temperature in degrees Celsius. Default is 75.
--max_temperature_difference : Maximal difference between minimum and maximum primer melting temperatures of selected primers. Default is 5.
--end_at_threshold     : Maximum allowed A/T nucleotides in final 3 nucleotides (3'-end). Default is 2.
--end_gc_threshold     : Maximum allowed G/C nucleotides in final 5 nucleotides (5'-end). Default is 3.
--monorun_threshold    : Maximum allowed length of a single nucleotide run. Default is 3.
--duorun_threshold     : Maximum allowed length of a double nucleotide run. Default is 3.
--mfe_threshold        : Minimum required MFE for determining hairpin formation risk. Default is -5.
--self_complementarity_threshold : Maximum primer-primer complementarity in worst alignment between primers. Default is 10.
##Greedy algorithm parameters
-amps                  : Number of amplicons to find. Default is 10.
##Sequence parameters
-n                     : Number of sequences to include. Default is -1 (all input sequences).
--min_characters       : Minimum number of characters in a sequence. Default is -1 (no minimum).
--max_degeneracy       : Maximum degeneracy of a sequence. Default is -1 (no maximum).
--max_n                : Maximum number of N characters in a sequence. Default is -1 (no maximum).
##System parameters
-c                     : Number of cores to use in multiprocessing mode. Default is 1 (no multiprocessing).
-sd                    : Random seed to use when selecting subset of input sequences (only if -n is provided and smaller than actual number of sequences). Default is 0.
```

If both a primer threshold file and primer property thresholds are given through the CLI, AmpliDiff will prioritize threshold given through the CLI.

AmpliDiff will output the following four files:
- logfile_X.txt : this file contains information from the greedy amplicon selection step of AmpliDiff (e.g. accepted and rejected amplicons)
- primers_X.fasta : this file contains all the selected primers in fasta format. Names are of the form AMPLICON_N_FY in case of a forward primer, or AMPLICON_N_RY in case of a reverse primer
- runtimes_X.txt : this file contains information on the runtimes of AmpliDiff steps, as well as an ordered list of the selected amplicons and corresponding primers
- sequences_included_X.txt : this file contains the identifiers of sequences that were considered
Here X denotes the seed that was used (even when using all genomes), N the rank of an amplicon (when it was selected, 1 for first), and Y an index for the primer (e.g. AMPLICON_2_F3 referes to the third forward primer for amplicon 2).

### Example
The ```example``` directory has a small example with some test data to run AmpliDiff. It also includes the expected output in order to check if your installation gives the correct results.
