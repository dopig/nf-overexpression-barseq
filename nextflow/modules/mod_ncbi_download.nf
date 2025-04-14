#!/usr/bin/env nextflow

/*
 * Given an assembly ID, downloads gzipped genomic GFF, FNA, & feature tables
 */

process downloadNcbiAssembly {
    tag { "Downloading assembly $assembly_id" }

    input:
        val assembly_id

    output:
        path "${assembly_id}_*_genomic.fna.gz", emit: fasta
        path "${assembly_id}_*_genomic.gff.gz", emit: gff
        path "${assembly_id}_*_feature_table.txt.gz", emit: feature_table
        path "species.txt", emit: species

    script:
    """
    /app/ncbi_download.sh "${assembly_id}"
    """
}
