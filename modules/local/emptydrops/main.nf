process EMPTYDROPS {
    tag   "${meta.id}"
    label 'process_medium'

    // Runs on the RAW matrix: the empty droplets are the ambient profile it
    // tests each barcode against. Emits per-barcode CALLS, not a matrix --
    // R cannot write .h5ad (docs/python-r-bridge.md).

    input:
    tuple val(meta), path(raw_h5ad)

    output:
    tuple val(meta), path("${meta.id}.ed_calls.csv"), emit: calls
    path "versions.yml"                             , emit: versions

    script:
    def lower  = params.emptydrops_lower  ?: 100
    def niters = params.emptydrops_niters ?: 10000
    def fdr    = params.emptydrops_fdr    ?: 0.01
    """
    emptydrops.R \\
        --raw ${raw_h5ad} \\
        --out ${meta.id}.ed_calls.csv \\
        --lower ${lower} --niters ${niters} --fdr ${fdr} --seed ${params.seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep="."))')
        DropletUtils: \$(Rscript -e 'cat(as.character(packageVersion("DropletUtils")))')
    END_VERSIONS
    """
}
