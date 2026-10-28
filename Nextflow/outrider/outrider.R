#' OUTRIDER autoencoder find q
#' Last update: 20-01-2025
#' Script create ods objects
#' Argument 1= input path ods
#' Argument 2= qfile
#' Argument 3= samplesheet
#' Argument 4= output path ods
#' Argument 5= output path res
#' Argument 6= gene annotation table (gene_annotation.tsv, made from the GTF)


library(OUTRIDER)
library(dplyr)
library(tidyr)
library(data.table)


if(.Platform$OS.type == "unix") {
    register(MulticoreParam(workers=min(1, multicoreWorkers())))
} else {
    register(SnowParam(workers=min(1, multicoreWorkers())))
}


args <- commandArgs(trailingOnly = TRUE)

# Open OUTRIDER dataset
ods <- readRDS(args[1])
qtable <- read.table(args[2], header=TRUE, sep="\t")
samplesheet <- fread(args[3])


# Get the encodingDimension with the highest evaluation Loss
opt_q <- qtable$encodingDimension[which.max(qtable$evaluationLoss)]

print(paste0("Optimal q is: ",opt_q))

rds_out_path <- args[4]
res_out_path <- args[5]
annotation_file <- args[6]

iter <- 15


ods <- controlForConfounders(ods, q=opt_q, iterations=iter)

ods <- computePvalues(ods, alternative="two.sided", method="BY")
ods <- computeZscores(ods)
res <- results(ods, all=TRUE)

# Output the outrider dataset file, and the results table. 
saveRDS(ods, file=rds_out_path)

# Filter data based on samples present in samplesheet
res <- res[res$sampleID %in% samplesheet$sampleID]
print(paste0("res: ",res ))

# Gene annotation from the same GTF that was used for counting (offline, no
# org.Hs.eg.db / TxDb versions involved).
anno <- fread(annotation_file)
idx <- match(res$geneID, anno$gene_id)
# fall back to the ID without version suffix (e.g. external counts from another release)
nover <- is.na(idx)
idx[nover] <- match(sub("\\..*$", "", res$geneID[nover]), anno$gene_id_nover)
message(sprintf("%.1f%% of the result genes were found in the annotation",
                100 * mean(!is.na(idx))))
if (all(is.na(idx))) {
    stop("None of the OUTRIDER genes is in the annotation table; ",
         "are the gene IDs from the same GTF?")
}

res$hgncSymbol <- anno$gene_name[idx]
res$chr   <- sub("^chr", "", anno$chr[idx])
res$start <- anno$start[idx]
res$end   <- anno$end[idx]

names(res)[names(res) == 'geneID'] <- 'EnsemblID'

# solve issue unimplemented type 'list' in 'EncodeElement'? https://stackoverflow.com/questions/24829027/unimplemented-type-list-when-trying-to-write-table
# as dataframe, because it complains about being an atomic vector??
res <- as.data.frame(apply(res,2,as.character))

write.table(res, paste("combined_samples", res_out_path, sep='_'), sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)

for (sampleid in unique(res$sampleID)){
    sample_out_path = paste(sampleid, res_out_path, sep='_')
    write.table(res[res$sampleID == sampleid,], sample_out_path, sep='\t', append = FALSE, row.names = FALSE, col.names = TRUE)
}
