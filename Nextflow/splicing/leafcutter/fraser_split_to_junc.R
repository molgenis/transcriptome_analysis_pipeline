#' FRASER split read counts -> LeafCutter .junc file (BED12, regtools layout)
#' Argument 1= <sampleID>.splitCounts.RDS (GRanges from countSplitReads)
#' Argument 2= output .junc file
#'
#' FRASER's introns are 1-based and inclusive, like STAR's SJ.out.tab, so the
#' same coordinate convention as in sj_to_junc.py is used (chromStart = start-2,
#' chromEnd = end+1, block sizes 1,1).
#' Unstranded counting gives strand "*"; LeafCutter drops "?" junctions, so "*"
#' is written as "+" (all junctions of unstranded data share one strand).

suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)

gr <- readRDS(args[1])
stopifnot(is(gr, "GRanges"), "count" %in% colnames(mcols(gr)))

strand_chr <- as.character(strand(gr))
strand_chr[strand_chr == "*"] <- "+"

start <- start(gr) - 2L
end   <- end(gr) + 1L
keep  <- start >= 0L & mcols(gr)$count > 0

junc <- data.frame(
    chrom       = as.character(seqnames(gr)),
    chromStart  = start,
    chromEnd    = end,
    name        = "JUNC",
    score       = as.integer(mcols(gr)$count),
    strand      = strand_chr,
    thickStart  = start,
    thickEnd    = end,
    itemRgb     = "255,0,0",
    blockCount  = 2L,
    blockSizes  = "1,1",
    blockStarts = paste0("0,", end - start - 1L),
    stringsAsFactors = FALSE
)[keep, ]

write.table(junc, args[2], sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE)
message(nrow(junc), " of ", length(gr), " junctions written to ", args[2])
