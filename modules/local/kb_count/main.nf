process KB_COUNT {
    tag   "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("kb_out")                                                   , emit: counts
    tuple val(meta), path("${meta.id}.metrics")                                      , emit: metrics
    path "versions.yml"                                                              , emit: versions

    script:
    def chem = meta.chemistry == 'SC3Pv2' ? '10xv2' : '10xv3'
    def args = task.ext.args ?: '--workflow standard --filter bustools'
    // kb-python bundles x86_64 binaries; on arm64 the container built kallisto
    // and bustools from source, so point kb at those.
    """
    kb count \\
        -i ${index}/index.idx \\
        -g ${index}/t2g.txt \\
        -x ${chem} \\
        -t ${task.cpus} \\
        -o kb_out \\
        --kallisto \$(command -v kallisto) \\
        --bustools \$(command -v bustools) \\
        ${args} \\
        ${reads[0]} ${reads[1]}

    # inspect.json carries the barcode-on-whitelist stats QC_ASSERT gates on;
    # run_info.json carries the pseudoalignment rate.
    mkdir -p ${meta.id}.metrics
    cp kb_out/inspect.json kb_out/run_info.json ${meta.id}.metrics/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kb-python: \$(kb --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        kallisto: \$(kallisto version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        bustools: \$(bustools version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
