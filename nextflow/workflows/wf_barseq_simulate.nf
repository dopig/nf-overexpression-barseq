#!/usr/bin/env nextflow

include { simulateReads } from '../modules/mod_barseq_simulate.nf'

params.samples_tsv = "$projectDir/../../shared/example/bobaseq_barseq_samples.tsv"
params.gff_path = "$projectDir/../results/ref/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
params.plasmid_json_path = "$projectDir/../results/library/plasmids.json"
params.random_seed = "None"
params.winner_count = 5
params.winner_strength = 5000
params.count_range = "0 1000"
params.plusminus = 1

workflow barseqSimulate {
    take:
    samples_tsv
    gff
    plasmid_json
    random_seed_param
    winner_count
    winner_strength
    count_range
    plusminus

    main:
    random_seed = random_seed_param != null ? random_seed_param : "None"
    simulateReads(samples_tsv, gff, plasmid_json, random_seed, winner_count, winner_strength, count_range, plusminus)

    emit:
    fastq = simulateReads.out.fastq
    winners = simulateReads.out.winners
}

workflow {
    if (!params.samples_tsv) {
        exit 1, "Error: The parameter 'samples_tsv' is required. Please specify it using '--samples_tsv <value>'."
    }
    barseqSimulate(file(params.samples_tsv), file(params.gff_path), file(params.plasmid_json_path), 
        params.random_seed, params.winner_count, params.winner_strength, params.count_range, params.plusminus)
}
