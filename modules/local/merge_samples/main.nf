process MERGE_SAMPLES {
    label 'process_high'

    input:
    path h5ads

    output:
    path "merged.h5ad"      , emit: h5ad
    path "merge_report.json", emit: report
    path "versions.yml"     , emit: versions

    script:
    """
    merge_samples.py --inputs ${h5ads} --out merged.h5ad --report merge_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$(python -c 'import anndata; print(anndata.__version__)')
    END_VERSIONS
    """
}
