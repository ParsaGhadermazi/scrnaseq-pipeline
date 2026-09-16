process QC_ASSERT {
    tag   "${meta.id}"
    label 'process_low'

    // Fails the run loudly on a bad quantification rather than letting a
    // garbage matrix flow downstream. Cross-checks the DECLARED chemistry
    // against the OBSERVED valid-barcode fraction -- a low fraction is the
    // signature of chemistry mis-declaration, which otherwise inflates counts
    // silently. See docs/lectures/01-fastq-to-counts.md section 9.

    input:
    tuple val(meta), path(metrics), path(readinfo)

    output:
    tuple val(meta), path("*.qc_flags.json"), emit: flags
    path "versions.yml"                     , emit: versions

    script:
    def strict = params.qc_strict == false ? '' : '--strict'
    """
    qc_assert.py \\
        --metrics ${metrics} \\
        --readinfo ${readinfo} \\
        --quantifier ${meta.quantifier} \\
        --sample ${meta.id} \\
        --min-valid-barcodes ${params.min_valid_barcode_fraction} \\
        --out ${meta.id}.qc_flags.json ${strict}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
