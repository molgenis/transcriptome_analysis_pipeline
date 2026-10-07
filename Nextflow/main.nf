/**
Nextflow main OUTRIDER, FRASER and MAE Workflow
author: T Niemeijer
**/
nextflow.enable.dsl=2

include { Outrider; OutriderCount; MergeOutridercounts; CreateOutriderDataset; OutriderOptim; MergeQfiles } from "./outrider/outrider"
include { Fraser; MergeCounts; FraserCount } from "./fraser/fraser"
include { MAEreadCounting; GetMAEresults; SplitVCFs; ConcatMAEResults } from "./MAE/MAE"
include { ResultsToHtml } from "./html_report/report"

workflow Outrider_Fraser_MAE_nf {
    /* 
    Gagneurlab Outrider Nextflow implementation 
    */

    // Start Channel from the samplesheet, getting the required info per sample for featureCounts.
    Channel
    .fromPath( params.samplesheet )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple( row.sampleID, row.bamFile, row.pairedEnd, row.strandSpecific ) }
    | OutriderCount
    | view { row }
    | collect
    | set { merge_ch }

    // Merge the separate counts from the merge channel and create a OUTRIDER Dataset
    merge_ch
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
    Gagneurlab Fraser Nextflow implementation with external counts.
    Since it was a bit tricky to implement this parallized like Outrider
    It makes sure that symbolic links to the *bam/bai files are included
    */
    bamfiles_ch = Channel
    .fromPath( params.samplesheet )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> row.bamFile }

    baifiles_ch = bamfiles_ch.map { bamFile -> "${bamFile}.bai" }

    FraserCount(params.fraser.frasercountsR, params.samplesheet, bamfiles_ch.collect(), baifiles_ch.collect())
    MergeCounts(params.extcounts.folder, params.fraser.mergescriptR, FraserCount.out, params.extcounts.amount_fraser)
    Fraser(params.samplesheet, MergeCounts.out, params.fraser.fraserR)

    readcount_ch = Channel
    .fromPath( params.samplesheet )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple( row.sampleID, row.vcf, row.bamFile ) }
    .filter { it[1] != null }
    .filter { it[1] != "NA" }

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
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t')
    .map { row -> row.genePanel }
    .filter { it && it != "NA" }
    .map { file(it, checkIfExists: true) }
    .unique() // In case file is used for multiple samples.

    // report
    genepanel_ch.mix(Outrider.out, Fraser.out, GetMAEresults.out)
    | collect
    | map {it -> tuple( params.samplesheet, it )}
    | ResultsToHtml

}


workflow Outrider_counts_nf {

    Channel
    .fromPath( params.samplesheet )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> row.expressionFile }
    | view { "expression file: ${it}" }
    | collect
    | set { merge_ch }

    // Merge the separate counts from the merge channel and create a OUTRIDER Dataset
    merge_ch
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
    
}


workflow {

    if (params.mode == 'counts') {

        Outrider_counts_nf()

    } else if (params.mode == 'bam') {

        Outrider_Fraser_MAE_nf()

    } else {

        error "Invalid mode parameter: ${params.mode}"
    }
}
