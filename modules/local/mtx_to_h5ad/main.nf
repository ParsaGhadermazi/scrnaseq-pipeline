process MTX_TO_H5AD {
    tag   "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(mtx_dir), path(source_h5ad), path(rho)

    output:
    tuple val(meta), path("${meta.id}.corrected.h5ad"), emit: h5ad
    path "versions.yml"                               , emit: versions

    script:
    """
    mtx_to_h5ad.py --mtx-dir ${mtx_dir} --source ${source_h5ad} --rho ${rho} \\
        --out ${meta.id}.corrected.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$(python -c 'import anndata; print(anndata.__version__)')
    END_VERSIONS
    """
}
