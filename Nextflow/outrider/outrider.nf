#!/usr/bin/env nextflow
/**
Nextflow OUTRIDER workflow
using Rsubreads featureCounts to get count matrices, uses Gagneur-lab OUTRIDER
for calling outliers. 
Results will be stored in a outrider dataset (*.rds) and a resultsfile (*.tsv)
**/
nextflow.enable.dsl=2

process OutriderCount {
    time '1h'
    memory '8 GB'
    cpus 4

    publishDir "$params.output/counts", mode: 'copy'

    input:
        tuple val(sampleID), path(bamFile), val(pairedEnd), val(strandSpecific)
    output:
        path "${sampleID}_outrider_counts.tsv"
    script:
        """
        ml "${params.RbundelBioconductor}"
        Rscript "${params.outrider.outridercountsR}" "${sampleID}" "${bamFile}" "${params.featurecounts.genes_gtf}" "${pairedEnd}" "${strandSpecific}"
        """
}

process MergeOutridercounts {
    time '30m'
    memory '16 GB'
    cpus 1
    
    publishDir "$params.output/counts", mode: 'copy'

    input:
        path inputFiles 
    output:
        path "merged_outrider_counts.txt"
    script:
        """
        ml "${params.RbundelBioconductor}"
        Rscript ${params.outrider.mergecountsR} ${inputFiles}
        """
}


process CreateOutriderDataset{
    time '30m'
    memory '12 GB'
    cpus 1

    publishDir "$params.output/outrider", mode: 'copy'

    input:
        tuple path(outriderCounts), path(samplesheet)
    output:
        tuple path("outrider.rds"), path("q_values.txt")

    script: 
        """
        ml "${params.RbundelBioconductor}"
        Rscript "${params.outrider.outriderDatasetR}" "${outriderCounts}" "${samplesheet}" "${params.extcounts.folder}" "${params.extcounts.amount_outrider}" "${params.featurecounts.genes_gtf}"
        """
}

process OutriderOptim{
    // Outrider optimize functions
    time '4h'
    memory '16 GB'
    cpus 1

    input:
        tuple path(outriderDataset), val(q_value)
    output:
        path "*.tsv" // file with encdims specific to this Q.

    script: 
        """
        ml "${params.RbundelBioconductor}"
        Rscript "${params.outrider.outriderOptimR}" "${outriderDataset}" "${q_value}"
        """
}

process MergeQfiles {
    // Outrider optimize functions
    time '10m'
    memory '1 GB'
    cpus 1
    
    publishDir "$params.output/outrider/optim", mode: 'copy'

    input:
        path inputFiles
    output:
        path "merged_q_files.tsv"
    script:
        """
        ml "${params.RbundelBioconductor}"
        Rscript ${params.outrider.mergeQFiles} ${inputFiles}
        """
}

process Outrider {
    time '4h'
    memory '32 GB'
    cpus 1

    publishDir "$params.output/outrider", mode: 'copy'

    input:
        tuple path(outriderDataset), path(qfile), path(samplesheet)
    output:
        path "*.rds"
        path "*.tsv"

    script: 
        """
        ml "${params.RbundelBioconductor}"
        Rscript ${params.outrider.outriderR} ${outriderDataset} ${qfile} ${samplesheet} final_outrider.rds result_table_outrider.tsv ${params.genomeBuild}
        """
}
