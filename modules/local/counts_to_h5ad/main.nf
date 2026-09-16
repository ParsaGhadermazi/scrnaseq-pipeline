process COUNTS_TO_H5AD {
    tag   "${meta.id}"
    label 'process_medium'

    // The Part 1 / Part 2 boundary. Whichever quantifier ran, everything
    // downstream sees one schema. Dispatch is on meta.quantifier.

    input:
    tuple val(meta), path(quant)

    output:
    tuple val(meta), path("${meta.id}.h5ad")    , emit: h5ad
    tuple val(meta), path("${meta.id}.raw.h5ad"), emit: raw, optional: true
    path "versions.yml"                         , emit: versions

    script:
    def obs = groovy.json.JsonOutput.toJson(
        meta.findAll { k, v -> !(k in ['id','quantifier','expect_cells']) && v != null } )
    def rounding = params.fractional_counts_policy ?: 'round'
    """
    counts_to_h5ad.py \\
        --quant ${quant} \\
        --quantifier ${meta.quantifier} \\
        --sample ${meta.id} \\
        --obs '${obs}' \\
        --rounding ${rounding} \\
        --out ${meta.id}.h5ad \\
        --out-raw ${meta.id}.raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        anndata: \$(python -c 'import anndata; print(anndata.__version__)')
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
