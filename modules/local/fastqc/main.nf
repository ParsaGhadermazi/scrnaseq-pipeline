process FASTQC {
    tag   "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(read2)      // R2 only -- see docs/lectures/01 section 11

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path "versions.yml"            , emit: versions

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
