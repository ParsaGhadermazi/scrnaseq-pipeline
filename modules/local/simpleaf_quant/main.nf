process SIMPLEAF_QUANT {
    tag   "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("af_quant")                  , emit: counts
    tuple val(meta), path("${meta.id}.metrics")        , emit: metrics
    path "versions.yml"                                , emit: versions

    script:
    // cr-like-em resolves multimapping UMIs by EM, which produces FRACTIONAL
    // counts. COUNTS_TO_H5AD applies an explicit, recorded rounding policy --
    // downstream tools that assume integers would otherwise misbehave silently.
    def args  = task.ext.args ?: '--resolution cr-like-em --expected-ori fw'
    def chem  = meta.chemistry == 'SC3Pv2' ? '10xv2' : '10xv3'
    """
    export ALEVIN_FRY_HOME=./af_home
    simpleaf set-paths

    # simpleaf can only INFER a t2g map from a splici index built with
    # --fasta + --gtf. Indexes built from a bare transcriptome carry no
    # inferable map, so locate one and pass it explicitly rather than failing
    # deep inside quant.
    T2G=""
    # params.index must be simpleaf's OUTPUT directory (the one holding both
    # index/ and ref/), not the piscem index dir itself. Nextflow stages only
    # the declared path, so a sibling ../ref would not exist in the work dir.
    for c in ref/t2g_3col.tsv ref/t2g.tsv t2g_3col.tsv t2g.tsv t2g.txt; do
        if [ -f "${index}/\$c" ]; then T2G="--t2g-map ${index}/\$c"; break; fi
    done
    if [ -z "\$T2G" ]; then
        echo "ERROR: no transcript-to-gene map found in or beside ${index}" >&2
        echo "  looked for: ref/t2g_3col.tsv, ref/t2g.tsv, t2g_3col.tsv, t2g.tsv" >&2
        echo "  --index must be simpleaf's output dir (containing index/ and ref/)" >&2
        exit 1
    fi

    simpleaf quant \\
        --reads1 ${reads[0]} --reads2 ${reads[1]} \\
        --index ${index}/index \\
        --chemistry ${chem} \\
        --unfiltered-pl \\
        --threads ${task.cpus} \\
        --output af_quant \\
        \$T2G ${args}

    mkdir -p ${meta.id}.metrics
    cp af_quant/simpleaf_quant_log.json ${meta.id}.metrics/ 2>/dev/null || true
    find af_quant -maxdepth 2 -name '*.json' -exec cp {} ${meta.id}.metrics/ \\; 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf --version | sed 's/simpleaf //')
        salmon: \$(salmon --version | sed 's/salmon //')
        alevin-fry: \$(alevin-fry --version | sed 's/alevin-fry //')
    END_VERSIONS
    """
}
