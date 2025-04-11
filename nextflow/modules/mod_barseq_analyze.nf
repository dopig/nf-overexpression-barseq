#!/usr/bin/env nextflow

/*
 * Analyze bobaseq barseq experiments
 */

process multiCodes {
    tag { "Counting barcodes in $set"}

    container 'multicodes'

    input:
        tuple val(id), val(set), val(output_prefix)
        path reads

    output:
        path "${output_prefix}.codes", emit: codes
        path "${output_prefix}.counts", emit: counts

    script:
    """
    /app/src/MultiCodes.pl -minQuality 0 -bs4 -index $id -out $output_prefix < $reads
    """
}

process bobaseqFitness {
    tag { "Analyzing experimental outcomes" }

    publishDir "results/barseq/fitness", mode: 'copy'
    container 'r-fitness'

    input:
        path samples_tsv
        path mapping_dir
        path feature_table
        path multicodes_collect

    output:
        path "fitness.Rimage", emit: r_image
        path "fitness.tsv", emit: tsv

    script:
    """
    /app/build_r_image.R \
        --samples $samples_tsv \
        --mapping_dir $mapping_dir \
        --feature_table_path $feature_table
    """
}


process graphFitness {
    tag { "Plotting $style outcomes" }

    publishDir "results/barseq/fitness", mode: 'copy'
    container 'r-plot-fitness'

    input:
        tuple val(style), path(tsv)
        path r_image

    output:
        path "$style", emit: images

    script:
    """
    /app/plot_fitness.R "$style" "$tsv" "$r_image"
    """
}
