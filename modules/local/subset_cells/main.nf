process SUBSET_CELLS {
    tag   "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(raw_h5ad), path(calls)

    output:
    tuple val(meta), path("${meta.id}.cells.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.id}.ed.json")   , emit: report
    path "versions.yml"                           , emit: versions

    script:
    """
    subset_cells.py --raw ${raw_h5ad} --calls ${calls} \\
        --out ${meta.id}.cells.h5ad --report ${meta.id}.ed.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$(python -c 'import anndata; print(anndata.__version__)')
    END_VERSIONS
    """
}
