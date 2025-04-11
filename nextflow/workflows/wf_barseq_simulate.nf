#!/usr/bin/env nextflow

include { simulateReads } from '../modules/mod_barseq_simulate.nf'

params.samples_tsv = "$projectDir/../../data/reference/bobaseq_barseq_samples.tsv"
params.gff_path = "$projectDir/../results/ref/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
params.plasmid_json_path = "$projectDir/../results/library/plasmids.json"

workflow barseqSimulate {
    take:
    samples_tsv
    gff
    plasmid_json

    main:
    simulateReads(samples_tsv, gff, plasmid_json)

    emit:
    fastq = simulateReads.out.fastq
    winners = simulateReads.out.winners
}

workflow {
    if (!params.samples_tsv) {
        exit 1, "Error: The parameter 'samples_tsv' is required. Please specify it using '--samples_tsv <value>'."
    }
    barseqSimulate(file(params.samples_tsv), file(params.gff_path), file(params.plasmid_json_path))
}
