process PRELIM_CLUSTER {
    tag   "${meta.id}"
    label 'process_medium'

    // THROWAWAY. Exists only because SoupX needs cluster labels to identify
    // genes a group should not express. Discarded immediately afterwards --
    // nothing downstream sees these labels.

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}.prelim_clusters.csv"), emit: clusters
    path "versions.yml"                                    , emit: versions

    script:
    """
    prelim_cluster.py --input ${h5ad} --out ${meta.id}.prelim_clusters.csv \\
        --resolution ${params.prelim_resolution} --n-pcs ${params.n_pcs} \\
        --n-hvg ${params.n_hvg} --seed ${params.seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
