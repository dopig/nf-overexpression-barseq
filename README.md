# OverEx-Library-Seq

Building a(n eventual nextflow) workflow for analyzing experiments using barseq on over expression plasmids.

## Environment Setup

### Micromamba

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
### Docker

You need to have Docker installed and running in order to simulate and map reads. To generate consensus reads after simulating PacBio sequencing, you will need to have run `./src/build_docker_image.sh ccs`. To map those reads with Boba-seq, you will need to have run `./src/build_docker_image.sh bobaseq`.

Note: BobaSeq uses usearch, and you may need to update `docker/bobaseq/Dockerfile` to match the currently available [usearch binaries](https://www.drive5.com/usearch/download.html).


## Running the script

You can run the actual script by adjusting some parameters in `./run_makelibrary.sh`.  The most important parameter to set up is the path to the reference/source genome fasta file.  The script expects you to have a .gff file in the same directory with the same file name.

## Processess

### Simulated Data for Library Mapping
First steps are to design a plasmid and generate simulated data from that plasmid, as though it were sequenced in ways that would contain both inserts and barcode for eventually mapping the two.
