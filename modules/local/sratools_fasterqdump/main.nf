process SRATOOLS_FASTERQDUMP {
    tag   "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), val(srr)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: '--split-files --include-technical'
    """
    # --include-technical is not optional: SRA flags the 10x barcode read as
    # technical, and without it fasterq-dump silently returns R2 only.
    prefetch --max-size u --progress ${srr}
    fasterq-dump ${args} --threads ${task.cpus} --outdir . ${srr}
    pigz -p ${task.cpus} *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //')
    END_VERSIONS
    """
}
