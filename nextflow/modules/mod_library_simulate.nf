#!/usr/bin/env nextflow

/*
 * Processes to simulate creating and PacBio sequencing a new bobaseq library
 */

process designLibrary {
    tag { "Designing simulated library ${name} from ${assembly.simpleName}" }

    publishDir 'results/library/simulation', mode: 'copy', pattern: '{log.txt, plasmids.json}'

    container 'py-simulate'

    input:
        path assembly
        val name
        val unique_barcode_count
        val library_size
        val coverage
        val random_seed

    output:
        path "pcrs.fasta", emit: pcrs
        path "plasmids.json", emit: plasmid_json
        path "log.txt", emit: log

    script:
    def randomSeedArg = random_seed != "None" ? "--random-seed ${random_seed}" : ""
    """
    python3 /app/simulate_library.py $assembly ${randomSeedArg}\
        --unique-barcodes $unique_barcode_count \
        --library-size $library_size \
        --coverage $coverage
    """
}

process simulatePacBioReads {
    tag "Simulating PacBio run for library mapping"

    container 'pbsim'

    input:
        path template
        val passes

    output:
        path "sequenced.bam", emit: bam

    script:
    """
    pbsim --strategy templ --template $template --pass-num $passes --prefix sequenced --method qshmm --qshmm /opt/QSHMM-RSII.model
    """
}
