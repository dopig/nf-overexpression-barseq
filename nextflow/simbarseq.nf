#!/usr/bin/env nextflow

// nextflow/main.nf

process simulate_reads {
    tag "simulate"

    container "pysimbarseq"

    input:
    path samples_tsv
    path results
    path indexes

    output:
    path "results/barseq/reads"

    script:
    """
    /app/simulate_barseq_reads.py --samples-tsv-path $samples_tsv --multiplex-index-tsv $indexes $results
    """
}

params.samples_tsv = null
params.results = "$projectDir/results"
params.indexes = "$projectDir/../shared/external/primers/barseq4.index2"

workflow {
    if (!params.samples_tsv) {
        exit 1, "Error: The parameter 'samples_tsv' is required. Please specify it using '--samples_tsv <value>'."
    }

    simulate_reads(file(params.samples_tsv), file(params.results), file(params.indexes))
}
