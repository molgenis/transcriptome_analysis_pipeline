#!/usr/bin/env nextflow
/**
Nextflow MAE workflow
Uses Gagneur-lab MAE derived implementation from the DROP pipeline.. 
Results will be stored in a resultsfile (*.tsv)
**/

nextflow.enable.dsl=2


process MAEreadCounting {
    time '4h'
    memory '4 GB'
    cpus 1
    
    input:
        tuple val(sampleID), path(vcf), path(bamFile), val(fasta), path(fastafolder)

    output:
        tuple val("${sampleID}"), path("*_maecounts.tsv")
    script:
    """
    ${CMD_MAE} bash -c '
    set -e
    set -u

    # From DROP:
    if samtools view -H "${bamFile}" | grep -q "@RG"; then
        printf "BAM contains read groups RG, continuing with ASEReadCounter...\n"
    else
        printf "%s\n" "ERROR: BAM file doesnt contain Read Group Tag RG"
        printf "RG doesnt exist, it can be added using:\n"
        printf "   gatk AddOrReplaceGroups -R /path/to/reference -I /your/input.bam -O /your/output.bam --QUIET true\n"
        printf " https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-\n"
        printf " Try rerunning this module using the BAM with RG tags\n"
        exit 1
    fi

    vcf="${vcf}"

    base_vcf="\${vcf%.vcf.gz}"

    contig="\${base_vcf##*.}"

    # Remove INFO field, it sometimes contains whitespace or invalid line lengths
    bcftools annotate -x "INFO" "\${vcf}" -o "noinfo_\${base_vcf}.vcf"

    # Create index
    gatk IndexFeatureFile -I "noinfo_\${base_vcf}.vcf"

    # Only select heterozygous variants and only SNPs
    gatk SelectVariants -V "noinfo_\${base_vcf}.vcf" -O "temp_\${base_vcf}.vcf" -R "${fasta}" --restrict-alleles-to BIALLELIC -select "vc.getHetCount()==1" --select-type-to-include SNP

    # Remove duplicates
    bcftools norm --rm-dup all "temp_\${base_vcf}.vcf" > "temp_\${base_vcf}.vcf.tmp" && mv "temp_\${base_vcf}.vcf.tmp" "temp_\${base_vcf}.vcf"

    # Create new index
    gatk IndexFeatureFile -I "temp_\${base_vcf}.vcf"

    # Run ASEReadCounter
    gatk ASEReadCounter -I "${bamFile}" -V "temp_\${base_vcf}.vcf" -L "\${contig}" -O "\${base_vcf}_temp_counts.tsv" -R "${fasta}"

    # Sorting the count file
    (head -n 1 "\${base_vcf}_temp_counts.tsv" && tail -n +2 "\${base_vcf}_temp_counts.tsv" | sort -k1,1 -V -s) > "\${base_vcf}_maecounts.tsv"
    '
    """
}

process GetMAEresults {

    time '10m'
    memory '4 GB'
    cpus 1

    publishDir "$params.output/MAE/DeSEQresults", mode: 'copy'

    input:
        val sampleid
        path asecounts
        path resultsR
        path annotateGeneR
        path annotation

    output:
        
        path "${sampleid}_result_mae.tsv"

    script:
        """
        ${CMD_MAE} Rscript "${resultsR}" "${asecounts}" "${sampleid}_result_mae_tmp.tsv" "${sampleid}"
        ${CMD_OUTRIDER_FRASER} Rscript "${annotateGeneR}" "${sampleid}_result_mae_tmp.tsv" "${sampleid}_result_mae.tsv" "${annotation}"
        """

}

process SplitVCFs {

    time '10m'
    memory '4 GB'
    cpus 1

    input:
        tuple val(sampleID), path(vcf), path(bamFile)

    output:
        tuple val(sampleID), path("split.${sampleID}.*.vcf.gz"), path(bamFile)


    script:
        """
        ${CMD_MAE} bash -c '
         # Index the VCF if not indexed already
        bcftools index ${vcf}

        # Get the list of chromosome names (skip the header)
        chroms=\$(bcftools view ${vcf} | grep -v '^#' | cut -f 1 | sort | uniq)

        for C in \${chroms}; do
            echo "Processing chromosome: \${C}"
            bcftools view -O z -o split.${sampleID}.\${C}.vcf.gz ${vcf} \${C}
        done
        '
        """
}


process ConcatMAEResults {

    time '5m'
    memory '4 GB'
    cpus 1

    publishDir "$params.output/MAE/ASEreadcounts", mode: 'copy'

    input:
        tuple val(sampleID), path(mae_files)

    output:
        val "${sampleID}"
        path "${sampleID}_result_mae.tsv"

    script:
    """
    files=( ${mae_files} )

    head -n 1 \${files[@]:0:1} > "tmp.${sampleID}_result_mae.tsv"

    for file in \${files[@]}; do
        tail -n+2 \$file >> "tmp.${sampleID}_result_mae.tsv"
    done

    cp "tmp.${sampleID}_result_mae.tsv" "${sampleID}_result_mae.tsv"

    """
}


