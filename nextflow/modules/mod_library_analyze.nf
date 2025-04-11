#!/usr/bin/env nextflow

/*
 * Get consensus from PacBio sequencing and map bobaseq library
 */

process runPbccs {
    tag { "Generating consensus of $name library reads" }

    publishDir 'results/library/analysis', mode: 'copy'

    container 'pbccs'

    input:
        path bam_path
        val name

    output:
        path "${name}.fastq.gz", emit: fastq // note: it's gzipped

    script:
    """
    ccs $bam_path "${name}.fastq"
    gzip "${name}.fastq"
    """
}

process bobaseqMap {
    tag { "Mapping $name library" }

    publishDir 'results/library/analysis', mode: 'copy'

    container 'bobaseq'

    input:
        val name
        path json
        path fasta
        path gff
        path oligos
        path fastq

    output:
        path "map_config.json", emit: json
        path "$name", emit: map
        path "Logs", emit: log

    script:

    """
    gunzip -c $fasta > /data/genome.fna
    gunzip -c $gff > /data/genome.gff
    gunzip -c "$fastq" > /data/fastq/${fastq.baseName}

    # Prepare modified config json file & then run bobaseq
    python3 /opt/bin/adapt_json.py $json --name $name --fasta genome.fna --gff genome.gff --oligos $oligos --out map_config.json
    python3 -W ignore::FutureWarning /app/src/run_steps.py map_config.json /data/fastq /work/map 1

    mv /work/map/* .
    """
}
