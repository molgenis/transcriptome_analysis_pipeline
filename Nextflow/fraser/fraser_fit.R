#' FRASER autoencoder fit + results (hyperparameter search done in fraser_optim_q.R)
#' Gagneur-lab FRASER (2.0)
#' Processes start from a samplesheet with SampleID's BAM paths featurecount settings etc. 
#' and creates fraser rds object and results .tsv
#' 28-10-2023 (updated 27-01-2025)
#' Argument 1: Samplesheet
#' Argument 2= Input/output 
#' Argument 3= Genome build (hg38 or hg19)

library(FRASER)
library(dplyr)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


args <- commandArgs(trailingOnly = TRUE)

# Setup parallelisation
if(.Platform$OS.type == "unix") {
    register(MulticoreParam(workers=min(4, multicoreWorkers())))
} else {
    register(SnowParam(workers=min(4, multicoreWorkers())))
}

workdir <- args[2]

genome_build <- args[3]

# Load original sample table

original_settingsTable <- fread(args[1])

fds <- loadFraserDataSet(dir=workdir)
fds <- filterExpressionAndVariability(fds, minDeltaPsi=0, filter=FALSE)
fds <- fds[mcols(fds, type="j")[,"passed"],]

# Hyperparameter search results: one <q>_hyperparams.tsv per tested q, made by
# fraser_optim_q.R (one job per q)
hp_files <- list.files(".", pattern = "_hyperparams\\.tsv$")
stopifnot(length(hp_files) > 0)
hp <- rbindlist(lapply(hp_files, fread), fill = TRUE)
hp <- hp[order(-aroc)]
write.table(hp, "fraser_hyperparameters.tsv", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
# every q of the grid must have been evaluated
if (file.exists("q_values.txt")) {
    not_done <- setdiff(as.numeric(readLines("q_values.txt")), hp$q)
    if (length(not_done) > 0) {
        stop("No hyperparameter results for q = ", paste(not_done, collapse = ", "))
    }
}
best_q <- hp$q[1]
message("Best q (highest aroc): ", best_q)
# keep the table in the dataset, so plotEncDimSearch() etc. keep working
metadata(fds)[["hyperParams_jaccard"]] <- hp
set.seed(42)
fds <- FRASER(fds, q=c(jaccard=best_q))


if (genome_build == "hg19") {
    fds <- annotateRanges(fds, GRCh = 37)
} else {
    fds <- annotateRanges(fds, GRCh = 38)
}


fds <- calculatePadjValues(fds, type="jaccard", geneLevel=TRUE) # geneLevel TRUE -> FALSE

saveFraserDataSet(fds, dir=workdir, name="fraser_out")

register(SerialParam())
res <- as.data.table(results(fds,aggregate=TRUE, all=TRUE))

res <- res[res$sampleID %in% original_settingsTable$sampleID] #filter out samplesheet samples

# Rename seqnames to chr
names(res)[names(res) == 'seqnames'] <- 'chr'

# Results per patient
for (sampleid in unique(res$sampleID)){
    sample_out_path = paste(sampleid, 'result_table_fraser.tsv', sep='_')
    write.table(res[res$sampleID == sampleid], sample_out_path, sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)
}

write.table(res, paste("combined_samples", 'result_table_fraser.tsv', sep='_'), sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)

