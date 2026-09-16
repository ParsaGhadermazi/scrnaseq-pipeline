process AMBIENT_SOUPX {
    tag   "${meta.id}"
    label 'process_medium'

    // The ONLY R step that returns a matrix, so the only one paying the MTX
    // conversion. Needs BOTH matrices: empties measure the soup.

    input:
    tuple val(meta), path(raw_h5ad), path(cells_h5ad), path(clusters)

    output:
    tuple val(meta), path("soupx_out")          , emit: mtx
    tuple val(meta), path("${meta.id}.rho.json"), emit: rho
    path "versions.yml"                         , emit: versions

    script:
    def force = params.soupx_rho ? "--rho-force ${params.soupx_rho}" : ''
    """
    ambient_soupx.R \\
        --raw ${raw_h5ad} --cells ${cells_h5ad} --clusters ${clusters} \\
        --outdir soupx_out --rho ${meta.id}.rho.json \\
        --min-rho ${params.soupx_min_rho} \\
        --tfidf-min ${params.soupx_tfidf_min} \\
        --soup-quantile ${params.soupx_soup_quantile} \\
        --seed ${params.seed} ${force}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SoupX: \$(Rscript -e 'cat(as.character(packageVersion("SoupX")))')
    END_VERSIONS
    """
}
