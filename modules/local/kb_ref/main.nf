process KB_REF {
    tag   "${prefix}"
    label 'process_high'

    // NOT VERIFIED against a real genome. KB_COUNT was tested with a
    // hand-built kallisto index; this module has not been run.
    storeDir "${params.index_cache}/kallisto"

    input:
    path fasta
    path gtf

    output:
    path "${prefix}"         , emit: index
    path "${prefix}/ref.json", emit: provenance
    path "versions.yml"      , emit: versions

    script:
    prefix = params.reference_id ?: "${fasta.baseName}_${gtf.baseName}_kb"
    // 'standard' gives spliced counts only. 'nac' (nascent/mature/ambiguous)
    // is kb's velocity-capable workflow and the analogue of a splici index.
    def wf = params.kb_workflow ?: 'standard'
    def extra = wf == 'nac'
        ? "-f2 ${prefix}/intron.fa -c1 ${prefix}/cdna_t2c.txt -c2 ${prefix}/intron_t2c.txt"
        : ""
    """
    mkdir -p ${prefix}

    # KB_COUNT expects index.idx and t2g.txt at the root of this directory.
    kb ref \\
        --workflow ${wf} \\
        -i ${prefix}/index.idx \\
        -g ${prefix}/t2g.txt \\
        -f1 ${prefix}/cdna.fa \\
        ${extra} \\
        --kallisto \$(command -v kallisto) \\
        --bustools \$(command -v bustools) \\
        ${fasta} ${gtf}

    cat <<-END_JSON > ${prefix}/ref.json
    {
      "reference_id": "${prefix}",
      "tool": "kb ref",
      "workflow": "${wf}",
      "genome_fasta": "${fasta}",
      "gtf": "${gtf}",
      "velocity_capable": ${wf == 'nac'}
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kb-python: \$(kb --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        kallisto: \$(kallisto version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
