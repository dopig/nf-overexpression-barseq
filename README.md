# OverEx-Library-Seq

Building a(n eventual nextflow) workflow for analyzing experiments using barseq on over expression plasmids.

## Environment Setup

# Micromamba

This project requires the following dependencies:

* Python 3.10
* NumPy
* SciPy
* Biopython
* samtools
* pbsim3

To create the environment, run the following commands:

```bash
micromamba create -f environment.yml -p ./envs/sim-overex-lib
micromamba activate ./envs/sim-overex-lib
```

# Docker

You need to have Docker installed and running in order to simulate reads.  To make the pbccs-container image, please run `./src/build_pbccs_image.sh` during set-up.

## Processess

### Simulated Data for Library Mapping
First steps are to design a plasmid and generate simulated data from that plasmid, as though it were sequenced in ways that would contain both inserts and barcode for eventually mapping the two.
