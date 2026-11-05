#' FRASER gather step 2: assemble the FraserDataSet and calculate PSI values
#' Gagneur-lab FRASER (2.0)
#' Argument 1= samplesheet (tsv, all samples)
#' Argument 2= split_ranges.RDS
#' Argument 3= number of CPUs
#' Argument 4= stranded counting (true/false)
#' Input : fraser_output/ (from the split merge step; a private copy),
#'         *.nonSplicedCounts.h5 (one per sample)
#' Output: fraser_output/ with the saved FraserDataSet (same layout as the
#'         old frasercounts.R produced)

suppressPackageStartupMessages({
    library(FRASER)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
rangesfile  <- args[2]
ncpu        <- as.integer(args[3])
stranded    <- tolower(args[4]) == "true"

register(MulticoreParam(workers = min(ncpu, multicoreWorkers())))

settings <- fread(samplesheet)
# only samples with a bam take part in FRASER (others may be OUTRIDER-only)
settings <- settings[!(is.na(settings$bamFile) | settings$bamFile %in% c("", "NA"))]
settings$bamFile <- basename(settings$bamFile)

# Strand protocol (0 = unstranded, 1 = forward, 2 = reverse) lives in a 'strand'
# column. For unstranded counting it is deliberately NOT set here: this dataset
# is later merged with unstranded external counts, whose samples have no strand
# value, and FRASER only accepts 0/1/2 in that column. A missing column means
# unstranded; FRASER's "Strand is not specified" warning is expected and harmless.
if (stranded) {
    settings$strand <- as.integer(settings$strandSpecific)
}
fds <- FraserDataSet(colData = settings, workingDir = "fraser_output")

# Place the per-sample results where FRASER expects its non-split cache.
h5files <- list.files(".", pattern = "\\.nonSplicedCounts\\.h5$")
ids <- sub("\\.nonSplicedCounts\\.h5$", "", h5files)
missing <- setdiff(samples(fds), ids)
if (length(missing) > 0) {
    stop("No non-split counts for sample(s): ", paste(missing, collapse = ", "))
}
for (i in seq_along(h5files)) {
    if (!ids[i] %in% samples(fds)) next
    stopifnot(file.copy(h5files[i],
                        FRASER:::getNonSplicedCountCacheFile(ids[i], fds),
                        overwrite = TRUE))
}

# Split counts: loaded from fraser_output/savedObjects (written by the merge
# step). Non-split counts: all samples hit the per-sample cache.
splitCounts <- getSplitReadCountsForAllSamples(fds = fds, recount = FALSE,
                                               BPPARAM = SerialParam())
splitCountRanges <- readRDS(rangesfile)
nonSplitCounts <- getNonSplitReadCountsForAllSamples(
    fds = fds, splitCountRanges = splitCountRanges, recount = FALSE,
    BPPARAM = SerialParam())

fds <- addCountsToFraserDataSet(fds, splitCounts, nonSplitCounts)

# Same as the tail of the old frasercounts.R
fds <- calculatePSIValues(fds)
fds <- saveFraserDataSet(fds)
