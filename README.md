# AmpliDiff

This repository contains the code for [AmpliDiff](https://www.biorxiv.org/content/10.1101/2023.07.22.550164v1), a Python tool that finds amplicons and corresponding primers in viral genomes in order to differentiate between different lineages/strains/species.

### Requirements
All dependencies can be found in the dependencies.txt files in the corresponding folders. Note that AmpliDiff has not been tested on Windows systems.

### Installation
AmpliDiff can simply be installed by cloning this repo, and building the amplicon_generation.pyx Cython file using the following commands:
```
git clone git@github.com:JaspervB-tud/AmpliDiff.git
cd AmpliDiff
python setup.py build_ext --inplace
```
As of now, AmpliDiff uses [Gurobi](https://www.gurobi.com) to solve the primer feasibility and minimization problems.
