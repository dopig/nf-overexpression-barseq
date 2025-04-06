#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process download_ncbi_assembly_files {
    tag { "Downloading assembly $assembly_id" }

    publishDir 'results/ref', mode: 'copy'

    container 'ncbi-image'

    input:
        val assembly_id

    output:
        path "${assembly_id}_*_genomic.fna.gz", emit: fasta
        path "${assembly_id}_*_genomic.gff.gz", emit: gff
        path "${assembly_id}_*_feature_table.txt.gz", emit: feature_table

    script:
    """
    docsum=\$(esearch -db assembly -query "$assembly_id" | efetch -format docsum)
    scaffold_count=\$(echo \$docsum | xtract -pattern Stat -if @category -equals scaffold_count -and @sequence_tag -equals all -element Stat)

    echo "Scaffold count: \${scaffold_count}"
    if [ \$scaffold_count != 1 ]; then
        echo "Script currently limited to organisms with 1 scaffold and $assembly_id has \$scaffold_count" >&2
        exit 1
    fi

    ftp_url=\$(echo \$docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)

    if [ -z "\$ftp_url" ]; then
        echo "Could not retrieve FTP URL for \$assembly_id" >&2
        exit 1
    fi

    # Extract the full file basename, e.g. GCF_000027325.1_ASM2732v1
    base_name=\$(basename "\$ftp_url")

    # Download the files
    for ext in genomic.fna.gz genomic.gff.gz feature_table.txt.gz; do
        wget -q "\${ftp_url}/\${base_name}_\${ext}"
    done
    """
}

process simulate_library {
    publishDir 'results/library', mode: 'copy'

    tag { "Simulating library ${name} from ${assembly.simpleName}" }

    container 'pysimlib'

    input:
        path pyscript_dir
        path assembly
        val name
        val unique_barcodes
        val library_size
        val coverage

    output:
        path "pcrs.fasta", emit: pcrs
        path "plasmids.json", emit: plasmids
        path "log.txt", emit: log

    script:
    """
    python3 $pyscript_dir/make_sim_library_chopped.py $assembly --unique-barcodes $unique_barcodes --library-size $library_size --coverage $coverage
    """
}


process simulate_pacbio_reads {
    container 'pbsim-image'

    input:
        path template
        path qshmm
        val passes

    output:
        path "sequenced.bam", emit: bam

    script:
    """
    pbsim --strategy templ --method qshmm --pass-num $passes --prefix sequenced --qshmm $qshmm --template $template
    """
}

process run_pbccs {
    publishDir 'results/library', mode: 'copy'

    container 'pbccs-image'

    input:
        path bam_path
        val name

    output:
        path "${name}.fastq", emit: fastq

    script:
    """
    ccs $bam_path "${name}.fastq"
    """
}

process bobaseq_map {
    publishDir 'results/map', mode: 'copy'
    container 'bobaseq'

    input:
        val name // library name
        path json
        path fasta
        path gff
        path oligos
        path fastq

    output:
        path "bobaseq_config.json", emit: json
        path "$name", emit: map
        path "Logs", emit: logs

    script:
    """
    echo "111"
    gunzip -c $fasta > genome.fna
    gunzip -c $gff > genome.gff

    echo "222"
    python3 /opt/bin/adapt_json.py $json --name $name --fasta genome.fna --gff genome.gff --oligos $oligos
    echo "333"
    mv ${fastq.name} /temp
    python3 /Boba-seq/src/run_steps.py bobaseq_config.json /temp /work/map 1
    mv /work/map/* .
    """
}


params.input = null
params.assembly_id="GCF_000027325.1"
params.library_name = "Mgenitalium"
params.oligos = "$projectDir/../../shared/reference/oligos.fasta"
params.unique_barcodes = 10000
params.library_size = 1000
params.coverage = 5

params.script_dir = "$projectDir/bin" // Temporary location for .py files that will eventually go into docker image/dockerfile
params.bobaseq_json="$projectDir/../../shared/reference/bobaseq_config.json"


workflow {
    download_ncbi_assembly_files(params.assembly_id)
    assembly = download_ncbi_assembly_files.out

    simulate_library(params.script_dir, assembly.fasta, params.library_name, params.unique_barcodes, params.library_size, params.coverage)

    simulate_pacbio_reads(
        simulate_library.out.pcrs,
        file('/Users/higgins/coding/bioinf-projects/pioneer/overex-library-seq/data/reference/pbsim-models/QSHMM-RSII.model'),
        15
    )

    run_pbccs(simulate_pacbio_reads.out.bam, params.library_name)

    bobaseq_map(params.library_name, params.bobaseq_json, assembly.fasta, assembly.gff, params.oligos, run_pbccs.out.fastq)
}
