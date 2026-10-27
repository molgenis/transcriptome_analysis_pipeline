#' FRASER autoencoder fit + results (hyperparameter search done in fraser_optim_q.R)
#' Gagneur-lab FRASER (2.0)
#' Processes start from a samplesheet with SampleID's BAM paths featurecount settings etc. 
#' and creates fraser rds object and results .tsv
#' 28-10-2023 (updated 27-01-2025)
#' Argument 1: Samplesheet
#' Argument 2= Input/output 
#' Argument 3= gene annotation table (gene_annotation.tsv, made from the GTF)
#' Argument 4= gene types to annotate with (comma separated, or 'all'), default protein_coding

library(FRASER)
library(dplyr)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)

# Setup parallelisation
if(.Platform$OS.type == "unix") {
    register(MulticoreParam(workers=min(4, multicoreWorkers())))
} else {
    register(SnowParam(workers=min(4, multicoreWorkers())))
}

workdir <- args[2]

annotation_file <- args[3]
biotypes <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else "protein_coding"

# ---- offline gene annotation ------------------------------------------------
gene_feature <- function(ranges, featureName, annotation){
    if (any(strand(ranges) == "*")) strand(annotation) <- "*"
    hits <- suppressWarnings(findOverlaps(ranges, annotation))
    dt <- data.table(from = queryHits(hits),
                     feature = mcols(annotation)[[featureName]][subjectHits(hits)])
    agg <- dt[, .(feature = paste(unique(feature), collapse = ";")), by = from]
    out <- rep(NA_character_, length(ranges))
    out[agg$from] <- agg$feature
    out
}

annotate_from_table <- function(fds, annofile, biotypes, featureName = "hgnc_symbol"){
    anno <- fread(annofile)
    if (!identical(biotypes, "all")) anno <- anno[gene_type %in% biotypes]
    anno <- anno[!is.na(gene_name) & gene_name != ""]
    stopifnot(nrow(anno) > 0)
    chr <- sub("^chr", "", anno$chr)
    if (all(startsWith(seqlevels(fds), "chr"))) chr <- paste0("chr", chr)
    gr <- GRanges(chr, IRanges(anno$start, anno$end), strand = anno$strand)
    mcols(gr)[[featureName]] <- anno$gene_name

    for (type in c("theta", "j")) {
        mcols(fds, type = type)[[featureName]] <-
            gene_feature(rowRanges(fds, type = type), featureName, gr)
    }
    share <- mean(!is.na(mcols(fds, type = "j")[[featureName]]))
    message(sprintf("%.1f%% of the junctions were annotated with a gene",
                    100 * share))
    if (share == 0) {
        stop("No junction could be annotated. Do the chromosome names of the ",
             "data and the GTF match, and does the GTF have gene_name?")
    }
    fds
}


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


fds <- annotate_from_table(fds, annotation_file, biotypes)


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

