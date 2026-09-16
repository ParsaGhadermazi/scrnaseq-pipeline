#!/usr/bin/env nextflow
/*
 * myRNASeqPipeline -- single-cell RNA-seq
 * Lecture 1: raw reads -> count matrices -> standardised .h5ad
 */
nextflow.enable.dsl = 2

include { SCRNASEQ } from './workflows/scrnaseq'

// ---------------------------------------------------------------- validation
def validateParams() {
    def errors = []

    if (!params.input) errors << "--input is required (samplesheet CSV)"
    if (!(params.stage in ['counts','qc','annotate','de']))
        errors << "--stage must be one of counts|qc|annotate|de, got '${params.stage}'"
    if (params.stage == 'qc' && !(params.qc_method in ['mad','fixed']))
        errors << "--qc_method must be 'mad' or 'fixed', got '${params.qc_method}'"
    if (!(params.step in ['fastq','counts']))
        errors << "--step must be 'fastq' or 'counts', got '${params.step}'"
    if (!(params.quantifier in ['cellranger','alevinfry','kallisto']))
        errors << "--quantifier must be one of cellranger|alevinfry|kallisto, got '${params.quantifier}'"

    if (params.stage == 'counts' && params.step == 'fastq') {
        if (params.quantifier == 'cellranger') {
            // cellranger is licensed and never shipped in the image
            if (!params.cellranger_path)
                errors << "--cellranger_path is required for quantifier 'cellranger': " +
                          "point at your unpacked 10x tarball (it is bind-mounted, never containerised)"
            else if (!file(params.cellranger_path).exists())
                errors << "--cellranger_path does not exist: ${params.cellranger_path}"
            if (!params.cellranger_ref && !(params.genome_fasta && params.gtf))
                errors << "quantifier 'cellranger' needs either --cellranger_ref (a prebuilt " +
                          "refdata-gex-* directory) or --genome_fasta and --gtf to build one with mkref"
        }
        else if (!params.index && !(params.genome_fasta && params.gtf)) {
            errors << "quantifier '${params.quantifier}' needs either --index, or --genome_fasta and --gtf to build one"
        }
    }

    if (errors) {
        error "Parameter validation failed:\n  - " + errors.join("\n  - ")
    }
}

def banner() {
    log.info """
    ================================================================
     ${workflow.manifest.name} v${workflow.manifest.version}
    ================================================================
     input        : ${params.input}
     outdir       : ${params.outdir}
     stage        : ${params.stage}
     step         : ${params.step}
     quantifier   : ${params.quantifier}
     chemistry    : ${params.chemistry ?: 'auto-detect from R1 length'}
     profile      : ${workflow.profile}
    ================================================================
    """.stripIndent()
}

// ---------------------------------------------------------------- entrypoint
workflow {
    validateParams()
    banner()
    SCRNASEQ()
}

// ---------------------------------------------------------- provenance record
// Nextflow only populates workflow.commitId when the pipeline is run AS A
// PROJECT (`nextflow run owner/repo`). Running a local `main.nf` -- which is how
// most people develop -- leaves it null, so the run manifest would record no
// code identity at all. Read git directly instead, so provenance works either
// way.
def readGitInfo() {
    def dir = workflow.projectDir.toString()
    def run = { cmd ->
        try {
            def p = ["bash", "-c", "git -C '${dir}' ${cmd} 2>/dev/null"].execute()
            p.waitFor()
            return p.exitValue() == 0 ? p.text.trim() : ''
        } catch (Exception e) { return '' }
    }
    def commit = run('rev-parse HEAD')
    if (!commit) {
        return [commit: 'n/a (not a git checkout)', branch: 'n/a',
                state: 'n/a', remote: 'n/a']
    }
    def dirty = run('status --porcelain')
    return [
        commit : commit,
        branch : run('rev-parse --abbrev-ref HEAD') ?: 'detached',
        // a dirty tree means the recorded commit does NOT describe what ran
        state  : dirty ? 'MODIFIED (uncommitted changes -- commit does not describe this run)' : 'clean',
        remote : run('config --get remote.origin.url') ?: 'none'
    ]
}
def gitInfo = readGitInfo()

workflow.onComplete {
    def dir = file("${params.outdir}/pipeline_info")
    dir.mkdirs()

    // Stage-specific provenance. Built here rather than inline: nesting a
    // triple-quoted string inside another is not valid Groovy.
    def analysisBlock = (params.stage == 'counts'
        ? [ "        stage            = ${params.stage}",
            "        step             = ${params.step}",
            "        quantifier       = ${params.quantifier}",
            "        chemistry        = ${params.chemistry ?: 'auto-detected per sample'}",
            "        reference        = ${params.cellranger_ref ?: params.index ?: params.genome_fasta ?: 'n/a'}" ]
        : [ "        stage            = ${params.stage}",
            "        counts_dir       = ${params.counts_dir ?: params.outdir + '/h5ad'}",
            "        emptydrops       = ${params.emptydrops_skip ? 'skipped' : 'lower=' + params.emptydrops_lower + ' niters=' + params.emptydrops_niters + ' fdr=' + params.emptydrops_fdr}",
            "        ambient_method   = ${params.ambient_method} (min_rho=${params.soupx_min_rho})",
            "        doublet_method   = ${params.doublet_method} (removal=${params.doublet_removal})",
            "        qc_method        = ${params.qc_method}" + (params.qc_method == 'mad' ? " nmads=${params.mad_nmads}" : " genes=${params.min_genes}-${params.max_genes} counts=${params.min_counts}-${params.max_counts}"),
            "        max_mito_pct     = ${params.max_mito_pct}",
            "        remove_hb_genes  = ${params.remove_hb_genes}",
            "        norm_method      = ${params.norm_method} (n_hvg=${params.n_hvg})",
            "        integration      = ${params.integration_method} (batch_key=${params.batch_key}, condition_key=${params.condition_key})" ]
        ).plus("        seed             = ${params.seed}").join('\n')

    file("${dir}/run_manifest.txt").text = """\
        # ${workflow.manifest.name} run manifest
        # Written ${new Date().format('yyyy-MM-dd HH:mm:ss')}

        [pipeline]
        version          = ${workflow.manifest.version}
        revision         = ${gitInfo.branch}
        commit           = ${gitInfo.commit}
        working_tree     = ${gitInfo.state}
        remote           = ${gitInfo.remote}
        project_dir      = ${workflow.projectDir}

        [execution]
        nextflow_version = ${nextflow.version}
        profile          = ${workflow.profile}
        container_engine = ${workflow.containerEngine ?: 'none'}
        container        = ${workflow.container ?: 'none'}
        command_line     = ${workflow.commandLine}
        launch_dir       = ${workflow.launchDir}
        work_dir         = ${workflow.workDir}
        start            = ${workflow.start}
        complete         = ${workflow.complete}
        duration         = ${workflow.duration}
        success          = ${workflow.success}
        exit_status      = ${workflow.exitStatus}

        [analysis]
${analysisBlock}
        [params]
        ${params.collect { k, v -> "${k} = ${v}" }.join('\n        ')}

        # Exact tool versions: pipeline_info/software_versions.yml
        # Exact environment : environment.lock.yml inside the container image
        """.stripIndent()

    log.info(workflow.success
        ? "\nDone. Manifest: ${dir}/run_manifest.txt\n"
        : "\nFailed after ${workflow.duration}. See ${dir}/run_manifest.txt\n")
}
