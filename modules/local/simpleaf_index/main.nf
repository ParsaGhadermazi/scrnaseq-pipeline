process SIMPLEAF_INDEX {
    tag   "${prefix}"
    label 'process_high'

    // NOT VERIFIED against a real genome. The module's interface was exercised
    // only via `simpleaf index --ref-seq` on a toy transcriptome; the splici
    // path below (--fasta + --gtf) has not been run.
    storeDir "${params.index_cache}/simpleaf"

    input:
    path fasta
    path gtf

    output:
    path "${prefix}"         , emit: index
    path "${prefix}/ref.json", emit: provenance
    path "versions.yml"      , emit: versions

    script:
    prefix = params.reference_id ?: "${fasta.baseName}_${gtf.baseName}_af"
    // spliced+intronic (splici) is the default deliberately. It keeps spliced
    // and unspliced counts separable, and that separation is the entire basis
    // of RNA velocity -- so this choice is what makes Part 13 (scVelo) possible
    // at all. A spliced-only index permanently forecloses it.
    def reftype = params.af_ref_type ?: 'spliced+intronic'
    def rlen    = params.read_length ?: 91
    """
    export ALEVIN_FRY_HOME=./af_home
    simpleaf set-paths

    # --rlen must match the biological read (R2) length: it sets how much
    # intronic flank is included around each exon-intron boundary.
    simpleaf index \\
        --output ${prefix} \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        --ref-type ${reftype} \\
        --rlen ${rlen} \\
        --threads ${task.cpus}

    cat <<-END_JSON > ${prefix}/ref.json
    {
      "reference_id": "${prefix}",
      "tool": "simpleaf index",
      "ref_type": "${reftype}",
      "read_length": ${rlen},
      "genome_fasta": "${fasta}",
      "gtf": "${gtf}",
      "velocity_capable": ${reftype.contains('intronic') || reftype.contains('unspliced')}
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf --version | sed 's/simpleaf //')
        piscem: \$(piscem --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
