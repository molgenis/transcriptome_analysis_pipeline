#' FRASER hyperparameter search, step 1: the grid of encoding dimensions (q)
#' Gagneur-lab FRASER (2.0)
#' Argument 1= FRASER dataset directory (read only)
#' Output: q_values.txt, one q per line (same log-spaced grid as FRASER's default)

suppressPackageStartupMessages(library(FRASER))

args <- commandArgs(trailingOnly = TRUE)

fds <- loadFraserDataSet(dir = args[1])
dontWriteHDF5(fds) <- TRUE

# Same filtering as in fraser_fit.R / the original fraser.R
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0, filter = FALSE)
fds <- fds[mcols(fds, type = "j")[, "passed"], ]

# FRASER's default range (getEncDimRange): 2 .. min(samples, junctions)/3,
# at most 15 log-spaced steps
a <- 2
b <- min(ncol(fds), nrow(fds)) / 3
if (b <= a) {
    q <- 2L
} else {
    q <- unique(round(exp(seq(log(a), log(b), length.out = min(15, b)))))
}
message("Encoding dimensions to test: ", paste(q, collapse = ", "))
writeLines(as.character(q), "q_values.txt")
