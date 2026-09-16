process DOUBLETS_SCDBLFINDER {
    tag   "${meta.id}"
    label 'process_medium'

    // Self-estimates the doublet rate from cell count, so unlike DoubletFinder
    // it needs neither nExp nor a pK sweep. Emits per-cell calls only.

    input:
    tuple val(meta), path(h5ad), path(clusters)

    output:
    tuple val(meta), path("${meta.id}.doublets_scdblfinder.csv"), emit: calls
    path "versions.yml"                                         , emit: versions

    script:
    def cl  = clusters.name != 'NO_FILE' ? "--clusters ${clusters}" : ''
    def dbr = params.doublet_rate ? "--dbr ${params.doublet_rate}" : ''
    """
    doublets_scdblfinder.R --input ${h5ad} \\
        --out ${meta.id}.doublets_scdblfinder.csv ${cl} ${dbr} --seed ${params.seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scDblFinder: \$(Rscript -e 'cat(as.character(packageVersion("scDblFinder")))')
    END_VERSIONS
    """
}
