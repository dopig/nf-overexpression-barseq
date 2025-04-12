#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { multiCodes } from '../modules/mod_barseq_analyze.nf'
include { bobaseqFitness } from '../modules/mod_barseq_analyze.nf'
include { graphFitness } from '../modules/mod_barseq_analyze.nf'

// Parameters
params.master_path = "$projectDir/../results"  // Default value (can be overridden with --master-path)

// Other parameters derived from master_path
params.reads = "${params.master_path}/barseq/reads/reads.fastq"
params.samples_tsv = "${params.master_path}/barseq/bobaseq_barseq_samples.tsv"
params.chosen_tsv = "${params.master_path}/barseq/reads/chosen_winners.tsv"

params.mapping_dir = "${params.master_path}/map/Mgenitalium"
params.feature_table = "${params.master_path}/ref/GCF_000027325.1_ASM2732v1_feature_table.txt.gz"


// You can opt to plot the chosen winners from a simulation rather than the top prots from fitness analysis
// If plot_all==true, simulation files will be grabbed from "${params.master_path}/barseq/reads/chosen_winners.tsv"
params.plot_all = false

workflow barseqAnalyze {
    take:
    samples_tsv
    fastq
    map
    feature_table
    chosen

    main:
    sample_ch = Channel.fromPath(samples_tsv)
        .splitCsv(header:true, sep:'\t')
        .map { row -> [row.index, row.set, row.sampleId + '_' + row.index]}
    multiCodes(sample_ch, fastq)

    bobaseqFitness(samples_tsv, map, feature_table, multiCodes.out.codes.collect())

    bobaseqFitness.out.tsv.view { f -> if (f.size() < 10) {"🟡 Bobaseq.R didn't find any top proteins."} }

    graph_input_ch = bobaseqFitness.out.tsv.filter { it.size() > 10 }.map { ['top-proteins', it] }
    if (chosen) {
        chosen_ch = chosen.map { ['chosen-winners', it] }
        graph_input_ch = graph_input_ch.mix(chosen_ch)
    }
    graphFitness(graph_input_ch, bobaseqFitness.out.r_image)
}

workflow {
    chosen = params.plot_all ? file(params.chosen_tsv) : null
    barseqAnalyze(file(params.samples_tsv), file(params.reads), file(params.mapping_dir), file(params.feature_table), chosen)
}
