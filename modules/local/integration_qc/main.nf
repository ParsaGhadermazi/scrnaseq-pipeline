process INTEGRATION_QC {
    label 'process_medium'

    input:
    path integrated

    output:
    path "integration_metrics.json", emit: metrics
    path "versions.yml"            , emit: versions

    script:
    """
    integration_qc.py --input ${integrated} --out integration_metrics.json \\
        --batch-key ${params.batch_key} --n-neighbors ${params.mixing_k} \\
        --subsample ${params.silhouette_subsample} --seed ${params.seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scikit-learn: \$(python -c 'import sklearn; print(sklearn.__version__)')
    END_VERSIONS
    """
}
