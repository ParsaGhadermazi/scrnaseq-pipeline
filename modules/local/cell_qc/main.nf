process CELL_QC {
    tag   "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(h5ad), path(doublet_calls)

    output:
    tuple val(meta), path("${meta.id}.qc.h5ad")     , emit: h5ad
    tuple val(meta), path("${meta.id}.qc_report.json"), emit: report
    path "versions.yml"                             , emit: versions

    script:
    def dbl = doublet_calls.name != 'NO_FILE' ? "--doublets ${doublet_calls}" : ''
    """
    cell_qc.py --input ${h5ad} ${dbl} \\
        --doublet-removal ${params.doublet_removal} \\
        --out ${meta.id}.qc.h5ad --report ${meta.id}.qc_report.json \\
        --method ${params.qc_method} --nmads ${params.mad_nmads} \\
        --min-genes ${params.min_genes} --max-genes ${params.max_genes} \\
        --min-counts ${params.min_counts} --max-counts ${params.max_counts} \\
        --max-mito-pct ${params.max_mito_pct} \\
        --min-cells-per-gene ${params.min_cells_per_gene} \\
        --remove-hb-genes ${params.remove_hb_genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
