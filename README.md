# nf-overexpression-barseq

**A modular Nextflow pipeline for simulating and analyzing barcoded overexpression libraries using barcode amplicon sequencing (barseq).**

The pipeline supports full workflows from synthetic library design to barcode mapping, experiment simulation, and fitness analysis, using containerized tools for portability and reproducibility.


## Project Structure

The pipeline is composed of four sub-workflows:

 1. Library Simulation: Design and barcoded plasmid library with random genomic inserts; simulate PacBio amplicon sequencing of the library with [pbsim3](https://github.com/yukiteruono/pbsim3).
 2.	Library Mapping: Get consensus sequences of library sequences with [ccs](https://ccs.how) and map barcodes to inserts with [Boba-seq](https://github.com/OGalOz/Boba-seq).
 3.	Barseq Simulation: Simulate experimental selection and barcode amplicon sequencing.
 4.	Barseq Analysis: Quantify barcodes with [MultiCodes](https://bitbucket.org/berkeleylab/feba) and compute relative fitness with [BobaseqFitness](https://github.com/morgannprice/BobaseqFitness).


## Requirements
 - [Nextflow (DSL2)](https://www.nextflow.io/docs/latest/install.html)
 - [Docker](https://www.docker.com/get-started/)

All dependencies are handled via Docker containers. No local installation of Python, Perl, or R packages required.


## Quick Start
1. Clone the repo
2. Run with 1 or 2 test parameters
    - This one should take about 30 seconds:
		```
		nextflow run main.nf -params-file shared/example/config/mu.yaml -profile local
		```
    - This one takes about 10 minutes:
		```
		nextflow run main.nf -params-file shared/example/config/pelagibacter-small.yaml -profile local
		```
3. Check out the outcomes in the results directory:
	```
	├── fitness
	│   ├── chosen-winners				How the actual simulated winner gene regions are plotted by Bobaseq
	│   ├── top_proteins.tsv			The top proteins in Bobaseq's opinion
	│   └── top-proteins				Bobaseq's plots of what it thinks are the top proteins
	├── library
	│   ├── logs						Logs from Bobaseq's analysis of the library
	│   └── plots						Plots from Bobaseq's analysis of the library
	├── nextflow_reports
	│   ├── report.html					Nextflow-generated performance metrics report
	│   └── trace.txt					Table summary of some of nextflow's performance metrics
	└── simulation
		├── chosen_winners.tsv			The actual simulated winner gene regions
		├── log-simulate-library.txt	Log from library simulation scripts
		└── log-simulate-reads.txt		Log from barseq read simulation scripts
	```


## Customization
Create a version of the yaml files to get the simulation that meets your needs by adjusting yaml and tsv files.

### YAML file
Here's how the example files (shared/examples/config/*.yaml) are set up:
```
# General
run_name: 			name for your results folder
library_name: 		name for the library
samples_tsv: 		path to spreadsheet summarizing experimental samples (described below)

# Simulation
random_seed:		setting this fixes the simulations from run to run

# Simulation - library
assembly_id:		ID of the NCBI assembly that will be downloaded for the simulation
unique_barcodes:	Size of the pool of barcodes available when building the library
library_size:		Size of the library constructed from the pool of barcodes
coverage: 5			Average number of reads per construct for PacBio library mapping
pbsim_passes:		Number of simulated passes PacBio makes over each DNA piece

# Simulation - barseq
winner_count:		Number of winning genes picked per condition
winner_strength:	How many more copies (on average) winner barcodes have over regular barcodes
count_range: 		String in format ("0 1000") for low/high range of read counts
plusminus:			How many reads replica samples can vary from each other
```
A few other parameters are standardized in nextflow.config, but you are much less likely to modify these.

### Samples file
These are samples for describing your experimental condition. These can be from real experiments! Example files are in shared/examples/samples/*.tsv. They are based on [Bt_bobaseq_exps_2023Sep11.tsv](https://figshare.com/ndownloader/files/42456846) for the [Bobaseq paper figshare page](https://figshare.com/articles/dataset/Barcoded_overexpression_screens_in_gut_Bacteroidales_identify_genes_with_new_roles_in_carbon_utilization_and_stress_resistance_/24195054?file=42456846).

The key columns are:
- index - The scripts expect this to match an index in [barseq4.index2](https://bitbucket.org/berkeleylab/feba/raw/0975b4c3239d35d041fa954fbc359e7f0cecea88/primers/barseq4.index2)
- t0set - This is used to group your experimental samples with appropriate t0 samples
- lib - Needs to have the same name as used in library_name in your yaml file
- Group & Condition - important for calling out T0 sample
- desc - Used by simulation scripts to group replica samples together


## Remote Use
The processes have been demonstrated to work well locally (Mac Silicon aka ARM64; use `-profile local`) and on an AWS EC2 *AMD64* instance (`-profile ec2`). [Ccs](https://ccs.how), which is used for getting the consensus of PacBio reads for library mapping, will not run on a ARM64 instance, though every other part of the pipeline will.
In addition to setting up AWS Batch (Compute Environment, etc.) and the appropriate IAM settings on AWS, you will need to make a few adaptations to run remotely on AWS.
- **AWS Security Credentials**: If you decide to run on AWS, you will need to modify `nextflow.config` to supply your specific information.
See [Nextflow: AWS Security Credentials](https://www.nextflow.io/docs/latest/aws.html#aws-security-credentials) for credential set-up. The file `nextflow.config` currently expects an AWS profile with the name 'nextflow-programmatic-access' to exist (option 2 from the link in the last sentence), but you can change or overpower that.
- **AWS Region**: The `us-west-2` region is currently hardcoded into nextflow.config, so you may want to change that.
- **S3 Bucket**: To run on AWS Batch from anywhere, you'll need to supply an S3 bucket for the Nextflow work directories. Additionally, if you run on EC2, `nextflow.config` is currently set up to save results to S3, so you'll want to supply the S3 bucket name, or change that setting if you run that way. You can either define your bucket in `nextflow.config` with `params.s3_bucket = "s3://YOUR_BUCKET"` or supply it via the `--s3_bucket` flag.

Note: the [Nextflow: AWS Custom AMI](https://www.nextflow.io/docs/latest/aws.html#custom-ami) guide will tell you to create a custom AMI (or install awscli in your Docker images) because the default AMI "Amazon ECS-Optimized Amazon Linux 2 (AL2) x86_64" lacks awscli, and I thought that was no longer important to worry about, because "Amazon ECS-Optimized Amazon Linux 2023 x86_64" already has awscli in it. That newer AMI works well for many things, but I was incapable of getting Nextflow to use the awscli it already had, so I needed to remove and reinstall it anyway (see: [AWS CLI docs](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)).


## Speeding Up the Workflow
The pbccs process is by far the most time consuming (>80% of the time spent) in my testing. Great time savings can be realized by running that process on a high CPU instance, and there is a flag for setting "--ccs_cpus" available for that purpose. Its memory needs increase with the CPUs allocated, so the memory requests are set to rise in tandem with CPU.


## Using with Real Data
Of course, analyzing _real_ data is even better! The named workflows `mapLibrary` and `barseqAnalyze` in this set up should significantly expedite analysis of real-world data.


## Find Out More
You can go to my website to find out more about this project and ways to contact me. Check out my post here: [Simulated Overexpression Screens: Building a Bioinformatics Pipeline Without Real Data](https://douglas-higgins.com/posts/bobaseq/).
