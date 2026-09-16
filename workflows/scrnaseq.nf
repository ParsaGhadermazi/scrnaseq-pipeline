/*
 * Top-level scRNA-seq workflow.
 * Lecture 1 scope. Parts 2+ attach after COUNTS_TO_H5AD.
 */

include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { SCRNASEQ_QC       } from '../subworkflows/local/scrnaseq_qc'
include { FASTQ_TO_COUNTS   } from '../subworkflows/local/fastq_to_counts'
include { COUNTS_TO_H5AD  } from '../modules/local/counts_to_h5ad/main'
include { QC_ASSERT       } from '../modules/local/qc_assert/main'
include { MULTIQC         } from '../modules/local/multiqc/main'

def parseRow(row) {
    // Columns the pipeline interprets. Anything else in the samplesheet is
    // carried verbatim into the meta map and lands in adata.obs -- so donor,
    // timepoint, tissue, sex all work with no schema change.
    // Declared INSIDE the method: a script-level `def` is local to run() and
    // is not visible from a method body.
    def reserved = ['sample','srr','fastq_1','fastq_2','chemistry','expect_cells','counts','batch']

    def meta = [:]
    meta.id           = row.sample
    meta.sample       = row.sample
    meta.chemistry    = row.chemistry?.trim()    ?: params.chemistry     // null => auto-detect
    meta.expect_cells = row.expect_cells?.trim() ?: params.expect_cells
    meta.batch        = row.batch?.trim()        ?: row.sample           // Part 3 uses this
    meta.quantifier   = params.quantifier

    // passthrough experimental design -> obs
    row.each { k, v -> if (!(k in reserved) && v?.trim()) meta[k] = v.trim() }
    return meta
}

workflow SCRNASEQ {

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    ch_rows = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    // ---- stage: qc -- start from per-sample .h5ad produced by stage 'counts' --
    // The samplesheet is the experiment manifest, reused across stages; files
    // are located by convention so you do not re-assemble paths between stages.
    if (params.stage == 'qc') {
        def dir = params.counts_dir ?: "${params.outdir}/h5ad"
        ch_qc_in = ch_rows.map { row ->
            def meta = parseRow(row)
            def f = file("${dir}/${meta.id}.h5ad")
            def r = file("${dir}/${meta.id}.raw.h5ad")
            if (!f.exists()) error "stage 'qc': missing ${f}\n  run --stage counts first, or set --counts_dir"
            if (!r.exists()) error "stage 'qc': missing ${r}\n  SoupX and emptyDrops need the RAW matrix; " +
                                   "it is not optional (see docs/python-r-bridge.md)"
            [ meta, f, r ]
        }
        SCRNASEQ_QC(ch_qc_in)
        ch_versions = ch_versions.mix(SCRNASEQ_QC.out.versions)

        ch_versions.unique().collectFile(name: 'software_versions.yml', sort: true,
                                         storeDir: "${params.outdir}/pipeline_info")
    }
    else {

    if (params.step == 'fastq') {

        // [ meta, srr ] or [ meta, [fastq_1, fastq_2] ]
        ch_input = ch_rows.map { row ->
            def meta = parseRow(row)
            if (row.srr?.trim()) {
                return [ meta, row.srr.trim(), [] ]
            }
            if (!row.fastq_1?.trim() || !row.fastq_2?.trim()) {
                error "Sample '${meta.id}': provide either 'srr' or both 'fastq_1' and 'fastq_2'"
            }
            return [ meta, null, [ file(row.fastq_1, checkIfExists: true),
                                   file(row.fastq_2, checkIfExists: true) ] ]
        }

        // Build the index if one was not supplied. Cached via storeDir.
        PREPARE_REFERENCE()
        ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

        FASTQ_TO_COUNTS(ch_input, PREPARE_REFERENCE.out.index)
        ch_counts    = FASTQ_TO_COUNTS.out.counts       // [ meta, filtered, raw ]
        ch_metrics   = FASTQ_TO_COUNTS.out.metrics      // [ meta, metrics.csv, readinfo.json ]
        ch_versions  = ch_versions.mix(FASTQ_TO_COUNTS.out.versions)
        ch_multiqc   = ch_multiqc.mix(FASTQ_TO_COUNTS.out.multiqc_files)

        QC_ASSERT(ch_metrics)
        ch_versions = ch_versions.mix(QC_ASSERT.out.versions)

    } else {
        // --step counts : start from existing matrices, skip all of Part 1
        ch_counts = ch_rows.map { row ->
            def meta = parseRow(row)
            if (!row.counts?.trim())
                error "Sample '${meta.id}': --step counts requires a 'counts' column " +
                      "pointing at the quantifier's output directory"
            [ meta, file(row.counts, checkIfExists: true) ]
        }
    }

    // Part 1 closes here: whichever quantifier ran, Part 2 sees one schema.
    COUNTS_TO_H5AD(ch_counts)
    ch_versions = ch_versions.mix(COUNTS_TO_H5AD.out.versions)

    if (!params.skip_multiqc) {
        MULTIQC(ch_multiqc.collect(sort: true).ifEmpty([]))
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    // one merged provenance file, deduplicated across parallel samples
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', sort: true,
                     storeDir: "${params.outdir}/pipeline_info")
    }

}
