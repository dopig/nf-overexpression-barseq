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
    echo "Horses"
    echo "$params.output_dir/$set"
    ${file(params.external_scripts_dir)}/src/MultiCodes.pl -minQuality 0 -bs4 -index $id -out $output_prefix < $params.reads
    """
}

process bobaseqFitness {
    publishDir "${file(params.output_dir)}", mode: 'copy'

    input:
        path multicodes_complete // Used only for dependency enforcement

    output:
        path "fitness.Rimage", emit: rimage
        path "fitness.tsv", emit: tsv

    script:
    """
    echo $params.samples_tsv
    echo $params.output_dir
    echo "NEver"

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
    publishDir "${file(params.output_dir)}/${exp_desc.trim().replaceAll(" +", "-")}", mode: 'copy'

    input:
        path r_image
        tuple val(exp_desc), val(locus_tag)

    output:
        path "${locus_tag}.svg", emit: image

    script:
    """
    #!/usr/bin/env Rscript

    load("${r_image}")
    svg("${locus_tag}.svg", width = 7, height = 7) # Adjust width and height in inches
    par(mar = c(4, 4.5, 2.5, 1), mgp = c(3, 1, 0), cex.main = 2, cex.axis = 1, cex.lab = 1.5)
    show("${exp_desc}", locus = "${locus_tag}", background = "Ec", ymax = 20,
        col.bg = "darkgrey", main = "${exp_desc}",
        extraLabels = data.frame(locus_tag = "${locus_tag}"))
    dev.off()
    """
}

// Parameters
params.external_scripts_dir = "../shared/external" // Location of external perl & R scripts
params.master_path = ''  // Default value (can be overridden with --master-path)

// Other parameters derived from master_path
params.reads = "${params.master_path}/barseq/reads/reads.fastq"
params.json_path = "${params.master_path}/barseq/reads/lib.json"
params.samples_tsv = "${params.master_path}/ref/bobaseq_barseq_samples.tsv"
params.output_dir = "${params.master_path}/barseq/fitness"

workflow {
    id_ch = Channel.fromPath(params.samples_tsv)
        .splitCsv(header:true, sep:'\t')
        .map { row -> [row.index, row.set, row.sampleId + '_' + row.index]}

    multiCodes(id_ch)

    bobaseqFitness(multiCodes.out.collect())

    graphFitness(
        bobaseqFitness.out.rimage,
        bobaseqFitness.out.tsv.splitCsv(sep:'\t', header:true)
                              .map{ [it.expDesc, it.locus_tag] }
    )
}
