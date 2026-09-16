/*
 * Stage `qc` -- Parts 2, 2-2 and 3.
 *
 * per-sample (parallel):
 *   EMPTYDROPS -> SUBSET_CELLS -> PRELIM_CLUSTER -> AMBIENT -> DOUBLETS -> CELL_QC
 * then merged:
 *   MERGE_SAMPLES -> NORMALIZE_HVG -> INTEGRATE -> INTEGRATION_QC
 *
 * The DAG has a LOOP in it: SoupX needs cluster labels, clustering needs
 * normalisation, but the real normalisation belongs after correction. Resolved
 * with a throwaway preliminary clustering (PRELIM_CLUSTER) whose labels are used
 * only by SoupX and never reach anything downstream.
 *
 * All joins are on meta.id, never on the whole meta map -- a map-valued join
 * silently yields an empty channel the moment any stage enriches meta.
 */

include { EMPTYDROPS           } from '../../modules/local/emptydrops/main'
include { SUBSET_CELLS         } from '../../modules/local/subset_cells/main'
include { PRELIM_CLUSTER       } from '../../modules/local/prelim_cluster/main'
include { AMBIENT_SOUPX        } from '../../modules/local/ambient_soupx/main'
include { MTX_TO_H5AD          } from '../../modules/local/mtx_to_h5ad/main'
include { DOUBLETS_SCDBLFINDER } from '../../modules/local/doublets_scdblfinder/main'
include { CELL_QC              } from '../../modules/local/cell_qc/main'
include { MERGE_SAMPLES        } from '../../modules/local/merge_samples/main'
include { NORMALIZE_HVG        } from '../../modules/local/normalize_hvg/main'
include { INTEGRATE            } from '../../modules/local/integrate/main'
include { INTEGRATION_QC       } from '../../modules/local/integration_qc/main'

workflow SCRNASEQ_QC {

    take:
    ch_input      // [ meta, filtered.h5ad, raw.h5ad ]

    main:
    ch_versions = Channel.empty()
    def NO_FILE = file("${projectDir}/assets/NO_FILE")

    ch_raw      = ch_input.map { meta, _f, raw -> [ meta, raw ] }
    ch_filtered = ch_input.map { meta, f, _r  -> [ meta, f   ] }

    // ---- 1. empty droplets -------------------------------------------------
    if (params.emptydrops_skip) {
        // Cell Ranger already applied an EmptyDrops-like call; re-running is optional.
        ch_cells = ch_filtered
    } else {
        EMPTYDROPS(ch_raw)
        ch_versions = ch_versions.mix(EMPTYDROPS.out.versions.first())

        SUBSET_CELLS(
            ch_raw.map { m, r -> [ m.id, m, r ] }
                  .join(EMPTYDROPS.out.calls.map { m, c -> [ m.id, c ] })
                  .map { _id, m, r, c -> [ m, r, c ] }
        )
        ch_versions = ch_versions.mix(SUBSET_CELLS.out.versions.first())
        ch_cells = SUBSET_CELLS.out.h5ad
    }

    // ---- 2. preliminary clustering (throwaway, for SoupX) ------------------
    PRELIM_CLUSTER(ch_cells)
    ch_versions  = ch_versions.mix(PRELIM_CLUSTER.out.versions.first())
    ch_clusters  = PRELIM_CLUSTER.out.clusters

    // ---- 3. ambient RNA correction ----------------------------------------
    if (params.ambient_method == 'none') {
        ch_corrected = ch_cells
    } else if (params.ambient_method == 'soupx') {
        AMBIENT_SOUPX(
            ch_raw.map { m, r -> [ m.id, m, r ] }
                  .join(ch_cells.map    { m, c -> [ m.id, c ] })
                  .join(ch_clusters.map { m, k -> [ m.id, k ] })
                  .map { _id, m, r, c, k -> [ m, r, c, k ] }
        )
        ch_versions = ch_versions.mix(AMBIENT_SOUPX.out.versions.first())

        MTX_TO_H5AD(
            AMBIENT_SOUPX.out.mtx.map { m, d -> [ m.id, m, d ] }
                .join(ch_cells.map            { m, c -> [ m.id, c ] })
                .join(AMBIENT_SOUPX.out.rho.map { m, j -> [ m.id, j ] })
                .map { _id, m, d, c, j -> [ m, d, c, j ] }
        )
        ch_versions  = ch_versions.mix(MTX_TO_H5AD.out.versions.first())
        ch_corrected = MTX_TO_H5AD.out.h5ad
    } else {
        error "ambient_method '${params.ambient_method}' is not implemented yet " +
              "(decontx and cellbender are planned -- see docs/COVERAGE.md)"
    }

    // ---- 4. doublets -------------------------------------------------------
    if (params.doublet_method == 'none') {
        ch_doublets = ch_corrected.map { m, _h -> [ m, NO_FILE ] }
    } else if (params.doublet_method == 'scdblfinder') {
        DOUBLETS_SCDBLFINDER(
            ch_corrected.map { m, h -> [ m.id, m, h ] }
                        .join(ch_clusters.map { m, k -> [ m.id, k ] })
                        .map { _id, m, h, k -> [ m, h, k ] }
        )
        ch_versions = ch_versions.mix(DOUBLETS_SCDBLFINDER.out.versions.first())
        ch_doublets = DOUBLETS_SCDBLFINDER.out.calls
    } else {
        error "doublet_method '${params.doublet_method}' is not implemented yet " +
              "(doubletfinder, scrublet, doubletdetection are planned)"
    }

    // ---- 5. cell + gene QC -------------------------------------------------
    CELL_QC(
        ch_corrected.map { m, h -> [ m.id, m, h ] }
                    .join(ch_doublets.map { m, d -> [ m.id, d ] })
                    .map { _id, m, h, d -> [ m, h, d ] }
    )
    ch_versions = ch_versions.mix(CELL_QC.out.versions.first())

    // ---- 6. merge + normalise + integrate ----------------------------------
    MERGE_SAMPLES(CELL_QC.out.h5ad.map { _m, h -> h }.collect(sort: true))
    NORMALIZE_HVG(MERGE_SAMPLES.out.h5ad)
    INTEGRATE(NORMALIZE_HVG.out.h5ad)
    INTEGRATION_QC(INTEGRATE.out.h5ad)
    ch_versions = ch_versions
        .mix(MERGE_SAMPLES.out.versions, NORMALIZE_HVG.out.versions,
             INTEGRATE.out.versions, INTEGRATION_QC.out.versions)

    emit:
    integrated  = INTEGRATE.out.h5ad          // the finalized merged table
    metrics     = INTEGRATION_QC.out.metrics
    per_sample  = CELL_QC.out.h5ad
    versions    = ch_versions
}
