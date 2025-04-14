#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { downloadNcbiAssembly } from '../modules/mod_ncbi_download.nf'
include { designLibrary } from '../modules/mod_library_simulate.nf'
include { simulatePacBioReads } from '../modules/mod_library_simulate.nf'

// Generally edited
params.assembly_id="GCF_000027325.1"
params.library_name = "Mgenitalium"
params.unique_barcodes = 10000
params.library_size = 500
params.coverage = 5
params.random_seed = "None"

// Sometimes edited
params.pbsim_passes = 15

workflow simulateLibrary {
    take:
    library_name
    assembly_id
    unique_barcodes
    library_size
    coverage
    pbsim_passes
    random_seed_param

    main:
    random_seed = random_seed_param != null ? random_seed_param : "None"
    downloadNcbiAssembly(assembly_id)
    designLibrary(downloadNcbiAssembly.out.fasta, library_name, unique_barcodes, library_size, coverage, random_seed)
    simulatePacBioReads(designLibrary.out.pcrs, pbsim_passes)

    emit:
    fasta = downloadNcbiAssembly.out.fasta
    gff = downloadNcbiAssembly.out.gff
    feature_table = downloadNcbiAssembly.out.feature_table
    plasmid_json = designLibrary.out.plasmid_json
    bam = simulatePacBioReads.out.bam
}

workflow {
    simulateLibrary(params.library_name, params.assembly_id, params.unique_barcodes, params.library_size,
        params.coverage, params.pbsim_passes, params.random_seed)
}
