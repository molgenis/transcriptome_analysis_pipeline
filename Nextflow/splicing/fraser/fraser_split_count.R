#' FRASER scatter step 1: split read counts for ONE sample
#' Gagneur-lab FRASER (2.0)
#' Argument 1= samplesheet (tsv, all samples)
#' Argument 2= sampleID to count
#' Argument 3= number of CPUs
#' Argument 4= stranded counting (true/false)
#' Output: <sampleID>.splitCounts.RDS (GRanges with split read counts)

suppressPackageStartupMessages({
    library(FRASER)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
samplesheet <- args[1]
sampleid    <- args[2]
ncpu        <- as.integer(args[3])
stranded    <- tolower(args[4]) == "true"

# Single-sample FraserDataSet. Same construction as the old frasercounts.R:
# the BAM is staged by Nextflow in the working directory, so use its basename.
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

# Counts the split reads and writes FRASER's own per-sample cache file.
# NcpuPerSample parallelises over chromosomes.
invisible(countSplitReads(sampleID = sampleid, fds = fds,
                          NcpuPerSample = ncpu, recount = TRUE))

# Hand the cache file to the next step under a predictable name.
cachefile <- FRASER:::getSplitCountCacheFile(sampleid, fds)
stopifnot(file.exists(cachefile))
stopifnot(file.copy(cachefile, paste0(sampleid, ".splitCounts.RDS")))
