process NORMALIZE_HVG {
    label 'process_high'

    input:
    path merged

    output:
    path "normalized.h5ad", emit: h5ad
    path "versions.yml"   , emit: versions

    script:
    """
    normalize_hvg.py --input ${merged} --out normalized.h5ad \\
        --method ${params.norm_method} --n-hvg ${params.n_hvg} \\
        --batch-key ${params.batch_key} --n-pcs ${params.n_pcs_merged} --seed ${params.seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
