/**
LeafCutter intron clustering + LeafcutterMD outlier splicing detection.

Junction input per sample, in order of preference:
  1. sjFile  (STAR SJ.out.tab)          -> SjToJunc
  2. bamFile (FRASER split read counts) -> SplitToJunc (made by FraserSplitCount)
Clustering is done once for all samples, the outlier test runs as one job per
chunk of clusters (the test of one cluster is independent of the others).
**/
nextflow.enable.dsl=2

process SjToJunc {
    time '10m'
    memory '2 GB'
    cpus 1

    input:
        tuple val(sampleID), path(sj)
        path script
    output:
        path "${sampleID}.junc"

    script:
        def multi = params.leafcutter.countMultiMapping ? '--with-multimapping' : ''
        """
        ${params.annotation.python} "${script}" "${sj}" --out "${sampleID}.junc" ${multi}
        """
}


process SplitToJunc {
    time '10m'
    memory { 4.GB * task.attempt }
    cpus 1

    input:
        path splitCounts      // <sampleID>.splitCounts.RDS
        path script
    output:
        path "*.junc"

    script:
        def sid = splitCounts.name.replace('.splitCounts.RDS', '')
        """
        ${CMD_OUTRIDER_FRASER} Rscript "${script}" "${splitCounts}" "${sid}.junc"
        """
}


process LeafcutterCluster {
    time { 4.h * task.attempt }
    memory { 16.GB * task.attempt }
    cpus 1

    publishDir "$params.output/leafcutterMD", mode: 'copy'

    input:
        path juncFiles
    output:
        path "leafcutter_perind_numers.counts.gz"

    script:
        """
        ls *.junc > juncfiles.txt
        ${CMD_LEAFCUTTER} python3 ${params.leafcutter.dir}/clustering/leafcutter_cluster_regtools.py \\
            -j juncfiles.txt -o leafcutter -r . \\
            -l ${params.leafcutter.maxIntronLen} -m ${params.leafcutter.minCluReads} -p ${params.leafcutter.minCluRatio}
        """
}


process LeafcutterSplit {
    time '10m'
    memory '4 GB'
    cpus 1

    input:
        path counts
        path script
    output:
        path "chunk_*.counts.gz"

    script:
        """
        ${params.annotation.python} "${script}" "${counts}" --chunks ${params.leafcutter.chunks} --prefix chunk
        """
}


process LeafcutterMD {
    time { 4.h * task.attempt }
    memory { 8.GB * task.attempt }
    cpus 2

    input:
        path chunk
    output:
        path "*_{pVals,effSize,clusterPvals}.txt"

    script:
        def prefix = chunk.name.replace('.counts.gz', '')
        """
        set +e
        ${CMD_LEAFCUTTER} Rscript ${params.leafcutter.dir}/scripts/leafcutterMD.R \\
            -o ${prefix} -p ${task.cpus} \\
            -s ${params.leafcutter.maxClusterSize} -c ${params.leafcutter.minCoverage} -t ${params.leafcutter.timeout} \\
            ${chunk} > md.log 2>&1
        status=\$?
        set -e
        cat md.log
        if [ \$status -ne 0 ]; then
            # leafcutterMD.R fails when no cluster of the chunk passes its filters; that is not an error here
            if grep -q "subscript out of bounds" md.log; then
                echo "No testable clusters in ${chunk}"
                for f in pVals effSize clusterPvals; do : > ${prefix}_\${f}.txt; done
            else
                exit \$status
            fi
        fi
        """
}


process LeafcutterResults {
    time '1h'
    memory { 8.GB * task.attempt }
    cpus 1

    publishDir "$params.output/leafcutterMD", mode: 'copy'

    input:
        path mdFiles
        path samplesheet
        path annotation
        path script
    output:
        path "*_result_table_leafcutterMD.tsv"
        path "leafcutterMD_*.txt.gz"

    script:
        """
        ${params.annotation.python} "${script}" --dir . --samplesheet "${samplesheet}" --annotation "${annotation}" --padj ${params.leafcutter.padj}
        """
}
