# transcriptome-analysis-pipeline 

This pipeline is based on OUTRIDER version 1.20.1 and FRASER version 1.99.4 (see "References").
These versions of OUTRIDER and FRASER were distributed under the MIT license. For more recent versions check out the corresponding repositories. 

## Getting started:

### 1. Clone git repository
```
git clone https://github.com/molgenis/transcriptome-analysis-pipeline.git
cd transcriptome-analysis-pipeline
```


### 2. Download images (on Sylabs cloud): 
```
# Add Sylabs Cloud and select it (one-time)
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud

# Recommended: point Apptainer cache to a fast local disk
export APPTAINER_CACHEDIR=/path/to/tmp

# Pull containers used by the pipeline
#   - FRASER / OUTRIDER (hg19/hg38 builds)
apptainer pull --dir "$APPTAINER_CACHEDIR" fraser_outrider.sif \
  library://timniem/rna_outliers/test:hg19hg38

#   - MAE
apptainer pull --dir "$APPTAINER_CACHEDIR" mae.sif \
  library://timniem/rna_outliers/test:tmae

#   - HTML report
apptainer pull --dir "$APPTAINER_CACHEDIR" report.sif \
  library://timniem/rna_outliers/test:report
```
Add these images to a folder called containers

### 3. Create a samplesheet (.tsv) with the following parameters:

| Column     | Description | Required | Default | 
| ----------- | ----------- | -------- | ------- | 
| sampleID    | Name of the sample | Yes | - | 
| bamFile   | path/to/the/bam (indexes should exist) | Yes  | - | 
| strandSpecific   | 0 (unstranded, 1 (forward stranded), 2 (reverse stranded) | Yes  | 0 | 
| pairedEnd   | TRUE or FALSE | Yes  | - | 
| excludeFit   | exclude this sample from autoencoder fit (TRUE or FALSE) | No  | FALSE |
| genePanel | hgnc genes in a .txt one gene per row, or a 3 column .bed file | No | - | 



## References

Brechtmann F, Mertes C, Matusevi??i??t?? A, et al. OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data. Am J Hum Genet. 2018;103(6):907-917. https://doi.org/10.1016/j.ajhg.2018.10.025

Scheller, I.F., Lutz, K., Mertes, C et al. Improved detection of aberrant splicing with FRASER 2.0 and the intron Jaccard index. Am Jrnl Hum Genet 110, 12 (2023). https://doi.org/10.1016/j.ajhg.2023.10.014

Yepez, V. A., Gusic, M., Kopajtich, R., Meitinger, T., Gagneur, J., & Prokisch, H. (2021). Gene expression and splicing counts from the Yepez, Gusic et al study - non strand-specific (1.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7126296

Y??pez, V. A., Smith, N. H., Mertes, C., & Gagneur, J. (2022). Gene expression and splicing counts from 49 tissues from GTEx v8 genome build hg38 - non-strand specific (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6078397
