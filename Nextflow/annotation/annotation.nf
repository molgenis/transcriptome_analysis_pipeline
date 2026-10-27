/**
Gene annotation table from the GTF that is used for counting.
Used for the OUTRIDER and the FRASER results (offline, no biomaRt / annotation packages).
**/
nextflow.enable.dsl=2

process GeneAnnotation {
    time '10m'
    memory '2 GB'
    cpus 1

    publishDir "$params.output/annotation", mode: 'copy'

    input:
        path gtf
        path script
    output:
        path "gene_annotation.tsv"

    script:
        """
        ${params.annotation.python} "${script}" "${gtf}" --out gene_annotation.tsv
        """
}
