#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process download_ncbi_assembly_files {
    tag { "Downloading assembly $assembly_id" }

    publishDir 'results/ref', mode: 'copy'

    container 'ncbi-download'

    input:
        val assembly_id

    output:
        path "${assembly_id}_*_genomic.fna.gz", emit: fasta
        path "${assembly_id}_*_genomic.gff.gz", emit: gff
        path "${assembly_id}_*_feature_table.txt.gz", emit: feature_table
        path "species.txt", emit: species

    script:
    """
    /app/ncbi_download.sh "${assembly_id}"
    """
}

process simulate_library {
    publishDir 'results/library', mode: 'copy'

    tag { "Simulating library ${name} from ${assembly.simpleName}" }

    container 'pysimlib'

    input:
        path pyscript_dir
        path assembly
        val name
        val unique_barcodes
        val library_size
        val coverage

    output:
        path "pcrs.fasta", emit: pcrs
        path "plasmids.json", emit: plasmids
        path "log.txt", emit: log

    script:
    """
    python3 /app/simulate_library.py $assembly --unique-barcodes $unique_barcodes --library-size $library_size --coverage $coverage
    """
}


process simulate_pacbio_reads {
    container 'pbsim'

    input:
        path template
        path qshmm
        val passes

    output:
        path "sequenced.bam", emit: bam

    script:
    """
    pbsim --strategy templ --method qshmm --pass-num $passes --prefix sequenced --qshmm $qshmm --template $template
    """
}

process run_pbccs {
    publishDir 'results/library', mode: 'copy'

    container 'pbccs'

    input:
        path bam_path
        val name

    output:
        path "${name}.fastq", emit: fastq

    script:
    """
    ccs $bam_path "${name}.fastq"
    """
}

process bobaseq_map {
    publishDir 'results/map', mode: 'copy'
    container 'bobaseq'

    input:
        val name // library name
        path json
        path fasta
        path gff
        path oligos
        path fastq

    output:
        path "map_config.json", emit: json
        path "$name", emit: map
        path "Logs", emit: logs

    script:

    """
    gunzip -c $fasta > /data/genome.fna
    gunzip -c $gff > /data/genome.gff
    mv ${fastq.name} /data/fastq

    # Prepare modified config json file & then run bobaseq
    python3 /opt/bin/adapt_json.py $json --name $name --fasta genome.fna --gff genome.gff --oligos $oligos --out map_config.json
    python3 -W ignore::FutureWarning /app/src/run_steps.py map_config.json /data/fastq /work/map 1

    mv /work/map/* .
    """
}

// Generally edited
params.assembly_id="GCF_000027325.1"
params.library_name = "Mgenitalium"
params.unique_barcodes = 500
params.library_size = 50
params.coverage = 10

// Sometimes edited
params.pbsim_passes = 15

// If you change this you need to change bobaseq_config.json (& make_sim_library.py if being uses) too
params.oligos = "$projectDir/../../shared/reference/oligos.fasta"

params.qshmm_path = "$projectDir/../../shared/reference/QSHMM-RSII.model"
params.script_dir = "$projectDir/bin" // Temporary location for .py files that will eventually go into docker image/dockerfile
params.bobaseq_json="$projectDir/../../shared/reference/bobaseq_config.json"

workflow {
    download_ncbi_assembly_files(params.assembly_id)
    assembly = download_ncbi_assembly_files.out

    simulate_library(params.script_dir, assembly.fasta, params.library_name, params.unique_barcodes, params.library_size, params.coverage)

    simulate_pacbio_reads(
        simulate_library.out.pcrs,
        file('/Users/higgins/coding/bioinf-projects/pioneer/overex-library-seq/data/reference/pbsim-models/QSHMM-RSII.model'),
        params.pbsim_passes
    )

    run_pbccs(simulate_pacbio_reads.out.bam, params.library_name)

    bobaseq_map(params.library_name, params.bobaseq_json, assembly.fasta, assembly.gff, params.oligos, run_pbccs.out.fastq)
}
