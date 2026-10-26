/**
Nextflow FRASER workflow

Counting is scattered per sample (like OutriderCount):
  FraserSplitCount      (per sample)  split reads
  FraserMergeSplit      (gather)      merge + intron filter
  FraserNonSplitCount   (per sample)  non-split reads at the filtered splice sites
  FraserFinalize        (gather)      build FraserDataSet, calculate PSI values
followed by MergeCounts (external counts) and Fraser (autoencoder + results).

Retry behaviour of the counting processes is set in nextflow.config (withName).
**/
nextflow.enable.dsl=2


process FraserSplitCount {
    time { 4.h * task.attempt }
    memory { 16.GB * task.attempt }
    cpus 2

    input:
        tuple val(sampleID), path(bamFile), path(baiFile)
        path samplesheet
        path splitcountR
    output:
        path "${sampleID}.splitCounts.RDS"

    script:
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${splitcountR}" "${samplesheet}" "${sampleID}" ${task.cpus} ${params.fraser.strandedCounting}
        """
}


process FraserMergeSplit {
    time { 4.h * task.attempt }
    memory { 64.GB * task.attempt }
    cpus 1

    input:
        path splitCounts      // all <sample>.splitCounts.RDS
        path bamFiles         // only the headers are read (chromosome naming style)
        path baiFiles
        path samplesheet
        path mergesplitR
        val minExpression
    output:
        path "fraser_output", emit: workdir
        path "split_ranges.RDS", emit: ranges

    script:
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${mergesplitR}" "${samplesheet}" ${minExpression} ${params.fraser.strandedCounting}
        """
}


process FraserNonSplitCount {
    time { 2.h * task.attempt }
    memory { 8.GB * task.attempt }
    cpus 4

    input:
        tuple val(sampleID), path(bamFile), path(baiFile)
        path ranges
        path samplesheet
        path nonsplitcountR
    output:
        path "${sampleID}.nonSplicedCounts.h5"

    script:
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${nonsplitcountR}" "${samplesheet}" "${sampleID}" "${ranges}" ${task.cpus} ${params.fraser.strandedCounting}
        """
}


process FraserFinalize {
    time { 4.h * task.attempt }
    memory { 48.GB * task.attempt }
    cpus 4

    input:
        path samplesheet
        path finalizeR
        path split_workdir, stageAs: "split_workdir"
        path ranges
        path nonSplitCounts   // all <sample>.nonSplicedCounts.h5
    output:
        path "fraser_output"

    script:
        """
        # work on a private copy: FRASER writes into its working directory
        cp -rL split_workdir fraser_output
        ${CMD_OUTRIDER_FRASER} Rscript "${finalizeR}" "${samplesheet}" "${ranges}" ${task.cpus} ${params.fraser.strandedCounting}
        """
}


// Legacy single-job counting (all samples in one task). No longer used by main.nf.
process FraserCount {
    time '10h'
    memory '32 GB'
    cpus 4

    input:
        path frasercountR
        path samplesheet
        path bamFiles
        path baiFiles
    output:
        path "fraser_output"

    script: 
        """
        mkdir "fraser_output"
        ${CMD_OUTRIDER_FRASER} Rscript "${frasercountR}" "${samplesheet}" "fraser_output"
        """
}


process MergeCounts {
    time '1h'
    memory '64 GB'
    cpus 1

    input:
        path ext_counts
        path mergescriptR
        path fraser_in, stageAs: "fraser_in"
        val ext_amount_fraser
        

    output:

        path "fraser_output"

    script: 
        """
        # work on a private copy: the R script saves into its input directory
        cp -rL fraser_in fraser_output
        ${CMD_OUTRIDER_FRASER} Rscript ${mergescriptR} "fraser_output" "${ext_counts}" "${ext_amount_fraser}"
        """

}

process FraserPrepareQ {
    // grid of encoding dimensions for the hyperparameter search
    time '2h'
    memory { 64.GB * task.attempt }
    cpus 1

    input:
        path fraser_in, stageAs: "fraser_in"
        path prepareqR
    output:
        path "q_values.txt"

    script:
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${prepareqR}" "fraser_in"
        """
}


process FraserOptimQ {
    // hyperparameter search for ONE encoding dimension (one job per q)
    time { 8.h * task.attempt }
    memory { 64.GB * task.attempt }
    cpus 2

    input:
        val q
        path fraser_in, stageAs: "fraser_in"
        path optimqR
    output:
        path "${q}_hyperparams.tsv"

    script:
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${optimqR}" "fraser_in" ${q}
        """
}


process Fraser {
    time '20h'
    memory '80 GB'
    cpus 6

    publishDir "$params.output/fraser", mode: 'copy'

    input:
        path samplesheet
        path fraser_in, stageAs: "fraser_in"
        path hyperparams      // all <q>_hyperparams.tsv
        path qvalues
        path fraserR

    output:

        path "*.tsv"

    script: 
        """
        # work on a private copy: the R script saves into its input directory
        cp -rL fraser_in fraser_work
        ${CMD_OUTRIDER_FRASER} Rscript "${fraserR}" "${samplesheet}" "fraser_work" "${params.genomeBuild}"
        """
}
