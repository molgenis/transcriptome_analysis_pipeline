/**
Nextflow main OUTRIDER, FRASER and MAE Workflow
author: T Niemeijer

One samplesheet drives everything. Per sample:
  - expressionFile given -> these OUTRIDER gene counts are used, no counting
  - no expressionFile    -> gene counts are made from the bamFile (featureCounts)
  - bamFile given        -> the sample is also used for FRASER (and MAE if a vcf is given)
  - no bamFile           -> the sample is only used for OUTRIDER
**/
nextflow.enable.dsl=2

include { Outrider; OutriderCount; OutriderQC; MergeOutridercounts; CreateOutriderDataset; OutriderOptim; MergeQfiles } from "./outrider/outrider"
include { Fraser; MergeCounts; FraserSplitCount; FraserMergeSplit; FraserNonSplitCount; FraserFinalize } from "./fraser/fraser"
include { MAEreadCounting; GetMAEresults; SplitVCFs; ConcatMAEResults } from "./MAE/MAE"
include { ResultsToHtml } from "./html_report/report"

// true when a samplesheet cell holds a value (not missing, empty or "NA")
def hasValue(v) {
    return v != null && v.toString().trim() != '' && v.toString().trim() != 'NA'
}

workflow Outrider_Fraser_MAE_nf {

    rows = file( params.samplesheet ).splitCsv( header: true, sep: '\t' )

    // Every sample needs counts from somewhere
    missing = rows.findAll { !hasValue(it.expressionFile) && !hasValue(it.bamFile) }
    if ( missing ) {
        error "These samples have neither an expressionFile nor a bamFile in the samplesheet: ${missing.collect { it.sampleID }.join(', ')}"
    }

    to_count   = rows.findAll { !hasValue(it.expressionFile) }
    precounted = rows.findAll {  hasValue(it.expressionFile) }
    fraser_rows = rows.findAll { hasValue(it.bamFile) }

    /*
    OUTRIDER: gene counts from the samplesheet when given, otherwise counted from the bam.
    */
    OutriderCount(
        Channel
        .fromList( to_count )
        .map { row -> tuple( row.sampleID, row.bamFile, row.pairedEnd, row.strandSpecific ) }
    )

    precounted_ch = Channel
    .fromList( precounted )
    .map { row -> file( row.expressionFile, checkIfExists: true ) }

    all_counts_ch = OutriderCount.out.counts.mix( precounted_ch )

    // Input QC: broken samples / wrong strandedness (see outrider/outrider_qc.py)
    OutriderQC(
        all_counts_ch.collect(),
        OutriderCount.out.summary.collect().ifEmpty( [] ),
        params.samplesheet,
        params.outrider.qcScript
    )

    // Merge the separate counts and create a OUTRIDER Dataset
    all_counts_ch
    .collect()
    | MergeOutridercounts
    | map { it -> tuple( it, params.samplesheet ) }
    | CreateOutriderDataset
    | set { optim_ch }

    // Find the optimal Q for the dataset.
    optim_ch
    .map { it -> it[1] }
    .splitCsv( header: false, sep: '\t' )
    .map { row -> tuple( "$params.output/outrider/outrider.rds", row )} //Hacky to include the outputdir outrider.rds.
    | OutriderOptim
    | collect
    | MergeQfiles
    | set { outrider_ch }

    // Start OUTRIDER
    outrider_ch
    | map { it -> tuple( "$params.output/outrider/outrider.rds", it, params.samplesheet )} //Hacky to include the outputdir outrider.rds.
    | Outrider

    /* 
    Gagneurlab Fraser Nextflow implementation with external counts, for the samples with a bam.
    Counting is scattered per sample, like OUTRIDER:
      split reads (per sample) -> merge + filter -> non-split reads (per sample)
      -> FraserDataSet + PSI values -> merge external counts -> fit
    */
    fraser_results = Channel.empty()

    if ( fraser_rows ) {
        if ( params.fraser.strandedCounting && params.extcounts.amount_fraser > 0 ) {
            error "fraser.strandedCounting=true cannot be combined with external FRASER counts (amount_fraser > 0): they are unstranded. Use strandedCounting=false or amount_fraser=0."
        }

        fraser_samples_ch = Channel
        .fromList( fraser_rows )
        .map { row -> tuple( row.sampleID, file(row.bamFile), file("${row.bamFile}.bai") ) }

        FraserSplitCount(fraser_samples_ch, params.samplesheet, params.fraser.splitCountR)

        FraserMergeSplit(
            FraserSplitCount.out.collect(),
            fraser_samples_ch.map { it[1] }.collect(),
            fraser_samples_ch.map { it[2] }.collect(),
            params.samplesheet,
            params.fraser.mergeSplitR,
            params.fraser.minExpressionInOneSample
        )

        FraserNonSplitCount(fraser_samples_ch, FraserMergeSplit.out.ranges.first(),
                            params.samplesheet, params.fraser.nonSplitCountR)

        FraserFinalize(params.samplesheet, params.fraser.finalizeR,
                       FraserMergeSplit.out.workdir, FraserMergeSplit.out.ranges,
                       FraserNonSplitCount.out.collect())

        MergeCounts(params.extcounts.folder, params.fraser.mergescriptR, FraserFinalize.out, params.extcounts.amount_fraser)
        Fraser(params.samplesheet, MergeCounts.out, params.fraser.fraserR)

        fraser_results = Fraser.out
    } else {
        log.warn "No sample has a bamFile: FRASER, MAE and the report are skipped (OUTRIDER only)."
    }

    /*
    MAE, for the samples with a vcf and a bam
    */
    readcount_ch = Channel
    .fromList( rows )
    .filter { hasValue(it.vcf) && hasValue(it.bamFile) }
    .map { row -> tuple( row.sampleID, row.vcf, row.bamFile ) }

    readcount_ch
    | SplitVCFs
    | set {splitvcf_ch}

    splitvcf_ch
    .flatMap { split_vcf ->
        split_vcf[1].collect { vcfFile ->
            tuple( split_vcf[0], vcfFile, split_vcf[2], params.fasta, params.fastafolder )
        }
    }
    | MAEreadCounting
    | set {merge_mae_ch}

    merge_mae_ch
    .groupTuple(by: 0)
    | ConcatMAEResults

    GetMAEresults(ConcatMAEResults.out, params.mae.resultsR, params.mae.annotateGeneR)

    genepanel_ch = Channel
    .fromList( rows )
    .map { row -> row.genePanel }
    .filter { hasValue(it) }
    .map { file(it, checkIfExists: true) }
    .unique() // In case file is used for multiple samples.

    // report (needs the FRASER results; samples without them are skipped in the report)
    if ( fraser_rows ) {
        genepanel_ch.mix(Outrider.out, fraser_results, GetMAEresults.out)
        | collect
        | map {it -> tuple( params.samplesheet, it )}
        | ResultsToHtml
    }

}


workflow {
    Outrider_Fraser_MAE_nf()
}
