/*
 * Part 1: reads -> count matrices.
 *
 * Quantifiers are swappable because each emits its NATIVE output directory as
 * [ meta, quant ]. COUNTS_TO_H5AD dispatches on meta.quantifier and is the only
 * place that knows about per-tool layouts. Adding STARsolo later is one module
 * file, one case here, and one reader function -- not a refactor.
 */

import groovy.json.JsonSlurper

include { SRATOOLS_FASTERQDUMP } from '../../modules/local/sratools_fasterqdump/main'
include { TENX_READ_DETECT     } from '../../modules/local/tenx_read_detect/main'
include { FASTQC               } from '../../modules/local/fastqc/main'
include { CELLRANGER_COUNT     } from '../../modules/local/cellranger_count/main'
include { SIMPLEAF_QUANT       } from '../../modules/local/simpleaf_quant/main'
include { KB_COUNT             } from '../../modules/local/kb_count/main'

workflow FASTQ_TO_COUNTS {

    take:
    ch_input      // [ meta, srr|null, [fastq_1, fastq_2]|[] ]
    ch_index      // reference index (provided, or built by PREPARE_REFERENCE)

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    ch_input
        .branch {
            from_sra  : it[1] != null
            from_local: true
        }
        .set { ch_branched }

    // --- fetch ---------------------------------------------------------------
    SRATOOLS_FASTERQDUMP(ch_branched.from_sra.map { meta, srr, _fq -> [ meta, srr ] })
    ch_versions = ch_versions.mix(SRATOOLS_FASTERQDUMP.out.versions.first())

    ch_raw_fastq = SRATOOLS_FASTERQDUMP.out.reads
        .mix(ch_branched.from_local.map { meta, _srr, fq -> [ meta, fq ] })

    // --- identify read roles + chemistry -------------------------------------
    TENX_READ_DETECT(ch_raw_fastq)
    ch_versions = ch_versions.mix(TENX_READ_DETECT.out.versions.first())

    // Fold what was DETECTED back into the meta map, so every downstream
    // process sees a concrete chemistry instead of a null. Parse once and
    // derive both channels, so reads and readinfo carry the SAME meta.
    ch_detected = TENX_READ_DETECT.out.reads
        .map { meta, reads, info ->
            def det = new JsonSlurper().parse(info.toFile())
            def files = (reads instanceof List ? reads : [reads]).sort { it.name }
            tuple( meta + [ chemistry: det.chemistry, r1_length: det.r1_length ], files, info )
        }
    ch_reads    = ch_detected.map { meta, files, _info -> [ meta, files ] }
    ch_readinfo = ch_detected.map { meta, _files, info -> [ meta, info  ] }

    // --- read QC : R2 only ---------------------------------------------------
    // FastQC on R1 does not crash, it lies: every read from one cell shares a
    // barcode, so duplication and overrepresentation go red by construction.
    if (!params.skip_fastqc) {
        FASTQC(ch_reads.map { meta, reads -> [ meta, reads[1] ] })
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(FASTQC.out.zip.map { _m, z -> z })
    }

    // --- quantify ------------------------------------------------------------
    switch (params.quantifier) {
        case 'cellranger':
            CELLRANGER_COUNT(ch_reads, ch_index)
            ch_counts   = CELLRANGER_COUNT.out.counts
            ch_metrics  = CELLRANGER_COUNT.out.metrics
            ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions.first())
            ch_multiqc  = ch_multiqc.mix(CELLRANGER_COUNT.out.metrics.map { _m, c -> c })
            break
        case 'alevinfry':
            SIMPLEAF_QUANT(ch_reads, ch_index)
            ch_counts   = SIMPLEAF_QUANT.out.counts
            ch_metrics  = SIMPLEAF_QUANT.out.metrics
            ch_versions = ch_versions.mix(SIMPLEAF_QUANT.out.versions.first())
            break
        case 'kallisto':
            KB_COUNT(ch_reads, ch_index)
            ch_counts   = KB_COUNT.out.counts
            ch_metrics  = KB_COUNT.out.metrics
            ch_versions = ch_versions.mix(KB_COUNT.out.versions.first())
            break
        default:
            error "Unknown quantifier '${params.quantifier}'"
    }

    // pair run metrics with detected read info so QC_ASSERT can cross-check the
    // declared chemistry against the observed valid-barcode fraction
    // Join on meta.id, not on the whole meta map: any stage that enriches meta
    // would otherwise silently break the join and skip QC entirely.
    ch_qc = ch_metrics.map    { meta, m -> [ meta.id, meta, m ] }
        .join( ch_readinfo.map { meta, i -> [ meta.id, i ] } )
        .map  { _id, meta, m, i -> [ meta, m, i ] }

    emit:
    counts        = ch_counts       // [ meta, quant ]  (native layout)
    metrics       = ch_qc           // [ meta, metrics, readinfo ]
    multiqc_files = ch_multiqc
    versions      = ch_versions
}
