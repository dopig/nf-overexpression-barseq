#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { runPbccs } from '../modules/mod_library_analyze.nf'
include { bobaseqMap } from '../modules/mod_library_analyze.nf'

params.library_name = null
params.bam = null
params.fasta = null
params.gff = null
params.bobaseq_json="$projectDir/../../shared/reference/bobaseq_config.json"

// If you change this you need to change bobaseq_config.json (& make_sim_library.py if being uses) too
params.oligos = "$projectDir/../../shared/reference/oligos.fasta"

workflow mapLibrary {
    take:
    library_name
    bam
    fasta
    gff
    bobaseq_json
    oligos

    main:
    runPbccs(bam, library_name)
    bobaseqMap(library_name, bobaseq_json, fasta, gff, oligos, runPbccs.out.fastq)

    emit:
    map = bobaseqMap.out.map
}

workflow {
    mapLibrary(params.library_name, params.bam, params.fasta,
        params.gff, params.bobaseq_json, params.oligos)
}
