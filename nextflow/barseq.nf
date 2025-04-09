#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process multiCodes {
    publishDir "$params.output_dir/$set", mode: 'copy'

    input:
        tuple val(id), val(set), val(output_prefix)

    output:
        path "${output_prefix}*", emit: files

    script:
    """
    ${file(params.external_scripts_dir)}/src/MultiCodes.pl -minQuality 0 -bs4 -index $id -out $output_prefix < $params.reads
    """
}

process bobaseqFitness {
    publishDir "${file(params.output_dir)}", mode: 'copy'

    input:
        path multicodes_complete // Used only for dependency enforcement

    output:
        path "fitness.Rimage", emit: r_image
        path "fitness.tsv", emit: tsv

    script:
    """
    build_r_image.R \
        --bobaseq_path ${file(params.external_scripts_dir)}/src/bobaseq.R \
        --samples ${file(params.samples_tsv)} \
        --json_path ${file(params.json_path)} \
        --codes_dir ${file(params.output_dir)} \
        --r_image fitness.Rimage \
        --output_tsv fitness.tsv
    """
}

process graphFitness {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(style), path(tsv)
        path r_image

    output:
        path "$style", emit: images

    script:
    """
    plot_fitness.R "$style" "$tsv" "$r_image"
    """
}


// Parameters
params.external_scripts_dir = "$projectDir/../shared/external" // Location of external perl & R scripts
params.master_path = "$projectDir/results"  // Default value (can be overridden with --master-path)

// Other parameters derived from master_path
params.reads = "${params.master_path}/barseq/reads/reads.fastq"
params.json_path = "${params.master_path}/barseq/reads/lib.json"
params.samples_tsv = "${params.master_path}/barseq/bobaseq_barseq_samples.tsv"
params.output_dir = "${params.master_path}/barseq/fitness"

// You can opt to plot the chosen winners from a simulation rather than the top prots from fitness analysis
// If plot_all==true, simulation files will be grabbed from "${params.master_path}/barseq/reads/chosen_winners.tsv"
params.plot_all = false

workflow {
    id_ch = Channel.fromPath(params.samples_tsv)
        .splitCsv(header:true, sep:'\t')
        .map { row -> [row.index, row.set, row.sampleId + '_' + row.index]}

    multiCodes(id_ch)

    bobaseqFitness(multiCodes.out.collect())

    graph_input_ch = bobaseqFitness.out.tsv.map { ['top-proteins', it] }
    if (params.plot_all) {
        chosen_tsv = "${params.master_path}/barseq/reads/chosen_winners.tsv"
        chosen_ch = Channel.fromPath(chosen_tsv).map { ['chosen-winners', it] }
        graph_input_ch = graph_input_ch.mix(chosen_ch)
    }

    graphFitness(graph_input_ch, bobaseqFitness.out.r_image)
}
