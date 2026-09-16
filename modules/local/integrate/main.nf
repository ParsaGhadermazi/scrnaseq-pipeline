process INTEGRATE {
    label 'process_high'

    // Corrects an EMBEDDING only. counts and X are untouched, so Part 5's DE
    // reads layers['counts'] and never sees corrected values.

    input:
    path normalized

    output:
    path "integrated.h5ad", emit: h5ad
    path "versions.yml"   , emit: versions

    script:
    """
    integrate.py --input ${normalized} --out integrated.h5ad \\
        --method ${params.integration_method} \\
        --batch-key ${params.batch_key} --condition-key ${params.condition_key} \\
        --n-pcs ${params.n_pcs} --seed ${params.seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')
        harmonypy: \$(python -c 'import harmonypy; print(harmonypy.__version__)' 2>/dev/null || echo 'n/a')
    END_VERSIONS
    """
}
