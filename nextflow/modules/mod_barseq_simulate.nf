#!/usr/bin/env nextflow

/*
 * Processes to simulate Illumina reads for a barseq experiment
 */

process simulateReads {
    tag "Simulating barseq experiment reads"

    input:
    path samples_tsv
    path gff
    path plasmid_json
    val random_seed
    val winner_count
    val winner_strength
    val count_range
    val plusminus

    output:
    path "chosen_winners.tsv", emit: winners
    path "log.txt", emit: log
    path "reads.fastq.gz", emit: fastq

    script:
    def randomSeedArg = random_seed != "None" ? "--random-seed ${random_seed}" : ""
    """
    python3 /app/simulate_barseq_reads.py ${randomSeedArg}\
        --samples-tsv-path $samples_tsv \
        --gff-path $gff \
        --plasmid-json-path $plasmid_json \
        --winner-count $winner_count \
        --winner-strength $winner_strength \
        --count-range $count_range \
        --plusminus $plusminus \
        --multiplex-index-tsv /app/barseq4.index2
    """
}
