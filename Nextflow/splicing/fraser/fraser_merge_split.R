#' FRASER gather step 1: merge per-sample split counts, filter introns
#' Gagneur-lab FRASER (2.0)
#' Argument 1= samplesheet (tsv, all samples)
#' Argument 2= minExpressionInOneSample (FRASER default 20)
#' Argument 3= stranded counting (true/false)
#' Input : *.splitCounts.RDS in the working directory (one per sample)
#' Output: fraser_output/ (merged split counts as HDF5) and split_ranges.RDS
#'         (introns that pass the filter; non-split reads are counted at these)

suppressPackageStartupMessages({
    library(FRASER)
    library(data.table)
    library(DelayedMatrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
minExpr     <- if (length(args) >= 2) as.numeric(args[2]) else 20
stranded    <- length(args) >= 3 && tolower(args[3]) == "true"

register(SerialParam())

settings <- fread(samplesheet)
# only samples with a bam take part in FRASER (others may be OUTRIDER-only)
settings <- settings[!(is.na(settings$bamFile) | settings$bamFile %in% c("", "NA"))]
settings$bamFile <- basename(settings$bamFile)

# FRASER takes the protocol from a 'strand' column (0 = unstranded, 1 = forward,
# 2 = reverse), not from 'strandSpecific'. Unstranded counting is required when
# unstranded external counts are merged in (params.fraser.strandedCounting = false).
settings$strand <- if (stranded) as.integer(settings$strandSpecific) else 0L
dir.create("fraser_output", showWarnings = FALSE)
fds <- FraserDataSet(colData = settings, workingDir = "fraser_output")

# Place the per-sample results where FRASER expects its split count cache.
rdsfiles <- list.files(".", pattern = "\\.splitCounts\\.RDS$")
ids <- sub("\\.splitCounts\\.RDS$", "", rdsfiles)
missing <- setdiff(samples(fds), ids)
if (length(missing) > 0) {
    stop("No split counts for sample(s): ", paste(missing, collapse = ", "))
}
for (i in seq_along(rdsfiles)) {
    if (!ids[i] %in% samples(fds)) next
    stopifnot(file.copy(rdsfiles[i],
                        FRASER:::getSplitCountCacheFile(ids[i], fds),
                        overwrite = TRUE))
}

# All samples hit the cache, so no BAM is read here (except for the header
# when mergeCounts checks the chromosome naming style).
splitCounts <- getSplitReadCountsForAllSamples(fds = fds, recount = FALSE,
                                               BPPARAM = SerialParam())

# Same filter as countRNAData(filter=TRUE)
maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpr
message(sum(passed), " of ", length(passed), " introns have >= ", minExpr,
        " split reads in at least one sample")
saveRDS(rowRanges(splitCounts[passed, ]), "split_ranges.RDS")
