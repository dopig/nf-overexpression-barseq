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

### Config Files

To run the [Boba-seq](https://github.com/OGalOz/Boba-seq) mapping, if your primers (around your barcode and around your insert) differ from the defaults, you will need to edit their sequences and/or names in `config/default_oligos.fasta` and give the script the locatoin where you save it. It's expecting it at `data/refs/oligos.fasta`. For reverse primers, you give the reverse complement of the physical primer sequence.
If these oligo names change (but not if the sequences change), you also need to change the names in `config/default_bobaseq_config.json`.

Your amplicon sequencing reaction (surrounding the insert and the barcode) will need a second set of primers.  If those, or their names, vary from defaults, you will also need to edit those in the `search_pcr2` section of `config/default_bobaseq_config.json`

Note: all off the values in all caps in `config/default_bobaseq_config.json` (e.g. "YOUR-GFF-FILENAME") will be changed programatically; do not worry about those.

## Running the script

You can run the actual script by adjusting some parameters in `./run_makelibrary.sh`.  The most important parameter to set up is the path to the reference/source genome fasta file.  The script expects you to have a .gff file in the same directory with the same file name.

## Processess

### Simulated Data for Library Mapping
First steps are to design a plasmid and generate simulated data from that plasmid, as though it were sequenced in ways that would contain both inserts and barcode for eventually mapping the two.
