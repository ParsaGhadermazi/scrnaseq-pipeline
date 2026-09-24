process CELLRANGER_COUNT {
    tag   "${meta.id}"
    label 'process_high'

    // cellranger is NOT in the image: 10x's licence forbids redistribution.
    // conf/base.config bind-mounts params.cellranger_path read-only at /opt/cellranger.
    // It is an x86_64 binary -- run with -profile amd on Apple silicon.

    input:
    tuple val(meta), path(reads)
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs")                               , emit: counts
    tuple val(meta), path("${meta.id}.metrics")                             , emit: metrics
    tuple val(meta), path("${meta.id}/outs/web_summary.html")              , emit: web_summary
    path "versions.yml"                                                    , emit: versions

    script:
    def args      = task.ext.args ?: '--create-bam=false'
    def chemistry = meta.chemistry ? "--chemistry=${meta.chemistry}" : '--chemistry=auto'
    def expect    = meta.expect_cells ? "--expect-cells=${meta.expect_cells}" : ''
    """
    export PATH=/opt/cellranger:\$PATH
    # Nextflow already staged these into the task dir; `cp -L` would duplicate
    # every FASTQ physically (~28 GB each for GSE174609) and burn inodes twice.
    # Relative symlinks resolve to the staged files, which are guaranteed
    # readable inside the container -- an absolute symlink could point outside
    # the bind-mounted paths and break under Apptainer.
    mkdir -p fastqs
    for f in ${reads}; do ln -sf "../\$f" "fastqs/\$f"; done

    cellranger count \\
        --id=${meta.id} \\
        --transcriptome=${reference} \\
        --fastqs=fastqs \\
        --sample=${meta.id} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${chemistry} ${expect} ${args}

    mkdir -p ${meta.id}.metrics
    cp ${meta.id}/outs/metrics_summary.csv ${meta.id}.metrics/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        reference: ${reference}
    END_VERSIONS
    """
}
