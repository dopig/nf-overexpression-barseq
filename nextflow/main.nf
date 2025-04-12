#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { simulateLibrary } from './workflows/wf_library_simulate.nf'
include { mapLibrary } from './workflows/wf_library_analyze.nf'
include { barseqSimulate } from './workflows/wf_barseq_simulate.nf'
include { barseqAnalyze } from './workflows/wf_barseq_analyze.nf'

//  Library simulation parameters - you surely want to change these
params.assembly_id="GCF_000027325.1"
params.library_name = "Mgenitalium"
params.unique_barcodes = 10000
params.library_size = 500
params.coverage = 5
params.pbsim_passes = 15
params.random_seed = 8 //"None"

// Barseq simulation parameters - you surely want to change this
params.samples_tsv = '../data/reference/bobaseq_barseq_samples.tsv'

// Library analysis parameter - some parameters can be tweaked
params.bobaseq_json="$projectDir/../shared/reference/bobaseq_config.json"
params.oligo_path = "$projectDir/../shared/reference/oligos.fasta"

workflow {
    simulateLibrary(params.library_name, params.assembly_id, params.unique_barcodes,
                    params.library_size, params.coverage, params.pbsim_passes, params.random_seed)

    mapLibrary(params.library_name, simulateLibrary.out.bam, simulateLibrary.out.fasta,
               simulateLibrary.out.gff, file(params.bobaseq_json), file(params.oligo_path))

    barseqSimulate(file(params.samples_tsv), simulateLibrary.out.gff,
                   simulateLibrary.out.plasmid_json, params.random_seed)

    barseqAnalyze(file(params.samples_tsv), barseqSimulate.out.fastq, mapLibrary.out.map,
                  simulateLibrary.out.feature_table, barseqSimulate.out.winners)
}
