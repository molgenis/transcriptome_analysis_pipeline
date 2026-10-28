#' FRASER scatter step 2: non-split read counts for ONE sample
#' Gagneur-lab FRASER (2.0)
#' Argument 1= samplesheet (tsv, all samples)
#' Argument 2= sampleID to count
#' Argument 3= split_ranges.RDS (filtered introns from the gather step)
#' Argument 4= number of CPUs
#' Argument 5= stranded counting (true/false)
#' Output: <sampleID>.nonSplicedCounts.h5

suppressPackageStartupMessages({
    library(FRASER)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
sampleid    <- args[2]
rangesfile  <- args[3]
ncpu        <- as.integer(args[4])
stranded    <- tolower(args[5]) == "true"

settings <- fread(samplesheet)
settings <- settings[settings$sampleID == sampleid]
stopifnot(nrow(settings) == 1)
settings$bamFile <- basename(settings$bamFile)

# FRASER takes the protocol from a 'strand' column (0 = unstranded, 1 = forward,
# 2 = reverse), not from 'strandSpecific'. Unstranded counting is required when
# unstranded external counts are merged in (params.fraser.strandedCounting = false).
settings$strand <- if (stranded) as.integer(settings$strandSpecific) else 0L

dir.create("fraser_output", showWarnings = FALSE)
fds <- FraserDataSet(colData = settings, workingDir = "fraser_output")

splitCountRanges <- readRDS(rangesfile)

# Rsubread::featureCounts at the donor/acceptor sites of the filtered introns;
# FRASER writes the result to its per-sample HDF5 cache.
invisible(countNonSplicedReads(sampleID = sampleid,
                               splitCountRanges = splitCountRanges,
                               fds = fds, NcpuPerSample = ncpu,
                               minAnchor = 5, recount = TRUE))

cachefile <- FRASER:::getNonSplicedCountCacheFile(sampleid, fds)
stopifnot(file.exists(cachefile))
stopifnot(file.copy(cachefile, paste0(sampleid, ".nonSplicedCounts.h5")))
