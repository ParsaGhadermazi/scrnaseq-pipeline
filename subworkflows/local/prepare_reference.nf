/*
 * Resolve the reference index for the selected quantifier.
 *
 * If an index was provided, use it. If not, BUILD it from genome_fasta + gtf.
 * Builds are cached with storeDir (see params.index_cache), so an index is
 * built once and reused across runs -- a human reference is hours of work and
 * many GB, and rebuilding it per run would make the pipeline unusable.
 */

include { CELLRANGER_MKREF } from '../../modules/local/cellranger_mkref/main'
include { SIMPLEAF_INDEX   } from '../../modules/local/simpleaf_index/main'
include { KB_REF           } from '../../modules/local/kb_ref/main'

workflow PREPARE_REFERENCE {

    main:
    ch_versions = Channel.empty()

    // cellranger keeps its own param because its reference is a 10x-specific
    // package, not an index in the same sense as the others.
    def provided = params.quantifier == 'cellranger' ? params.cellranger_ref : params.index

    if (provided) {
        ch_index      = Channel.value(file(provided, checkIfExists: true))
        ch_provenance = Channel.empty()
        log.info "Using provided reference: ${provided}"
    }
    else {
        ch_fasta = Channel.value(file(params.genome_fasta, checkIfExists: true))
        ch_gtf   = Channel.value(file(params.gtf,          checkIfExists: true))
        log.info "No index provided -- building one for '${params.quantifier}' " +
                 "(cached in ${params.index_cache})"

        switch (params.quantifier) {
            case 'cellranger':
                CELLRANGER_MKREF(ch_fasta, ch_gtf)
                ch_index      = CELLRANGER_MKREF.out.reference
                ch_provenance = CELLRANGER_MKREF.out.provenance
                ch_versions   = ch_versions.mix(CELLRANGER_MKREF.out.versions)
                break
            case 'alevinfry':
                SIMPLEAF_INDEX(ch_fasta, ch_gtf)
                ch_index      = SIMPLEAF_INDEX.out.index
                ch_provenance = SIMPLEAF_INDEX.out.provenance
                ch_versions   = ch_versions.mix(SIMPLEAF_INDEX.out.versions)
                break
            case 'kallisto':
                KB_REF(ch_fasta, ch_gtf)
                ch_index      = KB_REF.out.index
                ch_provenance = KB_REF.out.provenance
                ch_versions   = ch_versions.mix(KB_REF.out.versions)
                break
            default:
                error "Unknown quantifier '${params.quantifier}'"
        }
    }

    emit:
    index      = ch_index
    provenance = ch_provenance
    versions   = ch_versions
}
