/**
DESeq2 differential gene expression between the groups of the samplesheet
(column "group"), on the merged gene counts of the samples in the samplesheet.
**/
nextflow.enable.dsl=2

process Deseq2 {
    time '2h'
    memory { 16.GB * task.attempt }
    cpus 4

    publishDir "$params.output/deseq2", mode: 'copy'

    input:
        path counts
        path samplesheet
        path annotation
        path script
    output:
        path "deseq2_*"

    script:
        """
        ${CMD_MAE} Rscript "${script}" "${counts}" "${samplesheet}" "${annotation}" "${params.deseq2.reference}" "${params.deseq2.covariates}" ${params.deseq2.minCount} ${task.cpus}
        """
}
