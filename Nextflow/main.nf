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

include { Outrider; OutriderCount; OutriderQC; MergeOutridercounts; CreateOutriderDataset; OutriderOptim; MergeQfiles } from "./expression/outrider/outrider"
include { Fraser; FraserPrepareQ; FraserOptimQ; MergeCounts; FraserSplitCount; FraserMergeSplit; FraserNonSplitCount; FraserFinalize } from "./splicing/fraser/fraser"
include { MAEreadCounting; GetMAEresults; SplitVCFs; ConcatMAEResults } from "./MAE/MAE"
include { SjToJunc; SplitToJunc; LeafcutterCluster; LeafcutterSplit; LeafcutterMD; LeafcutterResults } from "./splicing/leafcutter/leafcutter"
include { Deseq2 } from "./expression/deseq2/deseq2"
include { GeneAnnotation } from "./annotation/annotation"
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

    // Gene annotation (one table from the counting GTF, for OUTRIDER and FRASER)
    GeneAnnotation( params.featurecounts.genes_gtf, params.annotation.script )

    // Start OUTRIDER
    outrider_final_ch = outrider_ch
    .map { it -> tuple( "$params.output/outrider/outrider.rds", it, params.samplesheet )} //Hacky to include the outputdir outrider.rds.

    Outrider( outrider_final_ch, GeneAnnotation.out )

    /* 
    Gagneurlab Fraser Nextflow implementation with external counts, for the samples with a bam.
    Counting is scattered per sample, like OUTRIDER:
      split reads (per sample) -> merge + filter -> non-split reads (per sample)
      -> FraserDataSet + PSI values -> merge external counts -> fit
    */
    fraser_results = Channel.empty()
    fraser_split_ch = Channel.empty()

    if ( fraser_rows ) {
        if ( params.fraser.strandedCounting && params.extcounts.amount_fraser > 0 ) {
            error "fraser.strandedCounting=true cannot be combined with external FRASER counts (amount_fraser > 0): they are unstranded. Use strandedCounting=false or amount_fraser=0."
        }

        fraser_samples_ch = Channel
        .fromList( fraser_rows )
        .map { row -> tuple( row.sampleID, file(row.bamFile), file("${row.bamFile}.bai") ) }

        FraserSplitCount(fraser_samples_ch, params.samplesheet, params.fraser.splitCountR)
        fraser_split_ch = FraserSplitCount.out

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

        // Hyperparameter search (encoding dimension q): one job per q, then the fit
        FraserPrepareQ(MergeCounts.out, params.fraser.prepareQR)

        q_ch = FraserPrepareQ.out
        .splitText()
        .map { it.trim() }
        .filter { it }

        FraserOptimQ(q_ch, MergeCounts.out.first(), params.fraser.optimQR)

        Fraser(params.samplesheet, MergeCounts.out, FraserOptimQ.out.collect(),
               FraserPrepareQ.out, params.fraser.fraserR, GeneAnnotation.out)

        fraser_results = Fraser.out
    } else {
        log.warn "No sample has a bamFile: FRASER, MAE and the report are skipped (OUTRIDER only)."
    }

    /*
    LeafcutterMD (outlier splicing). Junctions per sample come from
      - the STAR SJ.out.tab when the samplesheet has an sjFile, otherwise
      - the FRASER split read counts of the bam.
    */
    if ( params.leafcutter.enabled ) {
        lc_sj  = rows.findAll { hasValue(it.sjFile) }
        lc_bam = rows.findAll { !hasValue(it.sjFile) && hasValue(it.bamFile) }
        lc_n   = lc_sj.size() + lc_bam.size()

        if ( lc_n < params.leafcutter.minSamples ) {
            log.warn "LeafcutterMD skipped: ${lc_n} samples with junctions, at least ${params.leafcutter.minSamples} are needed (leafcutter.minSamples)."
        } else {
            SjToJunc(
                Channel.fromList( lc_sj ).map { row -> tuple( row.sampleID, file(row.sjFile, checkIfExists: true) ) },
                params.leafcutter.sjScript
            )

            // bam-only samples: FraserSplitCount has made <sampleID>.splitCounts.RDS for them
            lc_bam_ids = lc_bam.collect { it.sampleID + '.splitCounts.RDS' }
            SplitToJunc(
                fraser_split_ch.filter { it.name in lc_bam_ids },
                params.leafcutter.splitToJuncR
            )

            LeafcutterCluster( SjToJunc.out.mix( SplitToJunc.out ).collect() )
            LeafcutterSplit( LeafcutterCluster.out, params.leafcutter.splitScript )
            LeafcutterMD( LeafcutterSplit.out.flatten() )
            LeafcutterResults(
                LeafcutterMD.out.collect(),
                params.samplesheet,
                GeneAnnotation.out,
                params.leafcutter.resultsScript
            )
        }
    }

    /*
    DESeq2 differential expression between the groups of the samplesheet
    */
    if ( params.deseq2.enabled ) {
        groups = rows.findAll { hasValue(it.group) }.collect { it.group.toString().trim() }
        n_ok = groups.countBy { it }.findAll { it.value >= 2 }.size()
        if ( n_ok < 2 ) {
            log.warn "DESeq2 skipped: it needs a 'group' column with at least two groups of at least two samples."
        } else {
            Deseq2( MergeOutridercounts.out, params.samplesheet, GeneAnnotation.out, params.deseq2.script )
        }
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

    GetMAEresults(ConcatMAEResults.out, params.mae.resultsR, params.mae.annotateGeneR, GeneAnnotation.out)

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
