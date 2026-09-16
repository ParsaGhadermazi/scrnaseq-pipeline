process CELLRANGER_MKREF {
    tag   "${prefix}"
    label 'process_high'

    // NOT VERIFIED. cellranger mkref requires the licensed x86_64 binary and
    // far more RAM than a laptop has (~32 GB+, hours, ~11 GB output for human).
    // Written from 10x's documented interface; run it on a cluster first.
    //
    // storeDir, not publishDir: the index is an expensive, reusable artefact.
    // If ${prefix} already exists in the cache this process is skipped entirely.
    storeDir "${params.index_cache}/cellranger"

    input:
    path fasta
    path gtf

    output:
    path "${prefix}"        , emit: reference
    path "${prefix}.ref.json", emit: provenance
    path "versions.yml"     , emit: versions

    script:
    prefix = params.reference_id ?: "${fasta.baseName}_${gtf.baseName}_cr"
    // Biotype filtering is NOT cosmetic. cellranger counts a read only if it
    // maps uniquely to ONE gene; a read overlapping two annotated features is
    // discarded as ambiguous. Leaving pseudogenes in the GTF therefore does not
    // merely add junk rows -- it silently DELETES counts from the real genes
    // they resemble. See docs/lectures/01-fastq-to-counts.md section 7.
    def attrs = (params.gtf_biotypes as String).split(',')
                    .collect { "--attribute=gene_biotype:${it.trim()}" }.join(' \\\n        ')
    """
    export PATH=/opt/cellranger:\$PATH

    cellranger mkgtf ${gtf} filtered.gtf \\
        ${attrs}

    cellranger mkref \\
        --genome=${prefix} \\
        --fasta=${fasta} \\
        --genes=filtered.gtf \\
        --nthreads=${task.cpus} \\
        --memgb=${task.memory.toGiga()}

    # Two runs against different references are not comparable, so record what
    # this index actually is; COUNTS_TO_H5AD carries it into uns.
    cat <<-END_JSON > ${prefix}.ref.json
    {
      "reference_id": "${prefix}",
      "tool": "cellranger mkref",
      "genome_fasta": "${fasta}",
      "gtf": "${gtf}",
      "gtf_biotypes_kept": "${params.gtf_biotypes}",
      "cellranger_version": "\$(cellranger --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)"
    }
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
