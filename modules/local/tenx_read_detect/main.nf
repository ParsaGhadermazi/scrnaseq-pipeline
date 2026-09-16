process TENX_READ_DETECT {
    tag   "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(fastqs)

    output:
    // readinfo travels WITH the reads so the subworkflow can fold detected
    // chemistry into the meta map; also emitted alone for QC_ASSERT.
    tuple val(meta), path("renamed/*.fastq.gz"), path("*.readinfo.json"), emit: reads
    tuple val(meta), path("*.readinfo.json")                            , emit: readinfo
    path "versions.yml"                                                 , emit: versions

    script:
    def declared = meta.chemistry ? "--declared-chemistry ${meta.chemistry}" : ''
    """
    # Assign I1/R1/R2 by observed read length, not SRA's arbitrary _1/_2/_3
    # numbering, then rename to the 10x convention Cell Ranger requires.
    detect_10x_reads.py \\
        --fastqs ${fastqs} \\
        --sample ${meta.id} \\
        --outdir renamed \\
        --json ${meta.id}.readinfo.json \\
        ${declared}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        detect_10x_reads.py: \$(detect_10x_reads.py --version)
    END_VERSIONS
    """
}
