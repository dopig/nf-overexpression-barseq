#!/usr/bin/env nextflow

/*
 * Processes to simulate Illumina reads for a barseq experiment
 */

process simulateReads {
    tag "Simulating barseq experiment reads"

    publishDir 'results/barseq', mode: 'copy', pattern: '*.t*'

    container 'py-simbarseq'

    input:
    path samples_tsv
    path gff
    path plasmid_json
    val random_seed

    output:
    path "chosen_winners.tsv", emit: winners
    path "log.txt", emit: log
    path "reads.fastq", emit: fastq

    script:
    def randomSeedArg = random_seed != "None" ? "--random-seed ${random_seed}" : ""
    """
    /app/simulate_barseq_reads.py ${randomSeedArg}\
        --samples-tsv-path $samples_tsv \
        --gff-path $gff \
        --plasmid-json-path $plasmid_json \
        --multiplex-index-tsv /app/barseq4.index2
    """
}
