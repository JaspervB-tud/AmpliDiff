# Example

This folder contains some SARS-CoV-2 data from NCBIVirus as well as the results you should obtain by running AmpliDiff on it. The input genomes (pre-aligned) are found in `sequences_aligned.fasta`, and `metadata.tsv` contains the Pangolin lineage designations for these genomes.


To run AmpliDiff on this dataset, use the following command (assuming you are in this Example folder and you have the required dependencies installed):
```bash
python ../AmpliDiff/AmpliDiff.py sequences_aligned.fasta metadata.tsv -o Example_output/
```

Running the above command should give the following output:
```
Reading sequences
Done reading sequences
Randomly selecting up to -1 sequences with seed=0
Done selecting sequences
Processing sequences
Done processing sequences
Couldn't read primer_thresholds.csv, using default thresholds or those supplied as command line arguments
Generating primer index
Initially contains 39096 forward primers
Finally contains 3978 forward primers
Removed 0 primers occurring both as forward and reverse
Initially contains 39096 reverse primers
Finally contains 4012 reverse primers
Removed 0 primers occurring both as forward and reverse
Done generating primer index
Determining amplicon differentiabilities
Transforming input sequences to numeric representations
Done transforming input sequences
Determining sequence pairs with different lineages
Done determining pairs
Calculating amplicon differentiabilities
Done calculating differentiabilities
Done determining amplicon differentiabilities
Running greedy algorithm
...
Done running greedy algorithm
```
with a list of primers before the final line. The primers can now be found in `Example_output/primers_0.fasta`. As can be seen, there are only 4 amplicons (and a total of 16 primers) since after the 4th amplicon every pair of sequences with different lineages can be discriminated.
