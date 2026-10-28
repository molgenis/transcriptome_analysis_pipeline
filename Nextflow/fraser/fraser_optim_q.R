#' FRASER hyperparameter search, step 2: evaluate ONE encoding dimension (q)
#' Gagneur-lab FRASER (2.0)
#' Argument 1= FRASER dataset directory (read only)
#' Argument 2= q to evaluate
#' Output: <q>_hyperparams.tsv (q, noise, nsubset, eval, aroc)
#'
#' Every job uses the same seed before the outlier injection, so all q values
#' are evaluated on the same injected outliers and the results are comparable.

suppressPackageStartupMessages(library(FRASER))

args <- commandArgs(trailingOnly = TRUE)
q <- as.numeric(args[2])

register(SerialParam())

fds <- loadFraserDataSet(dir = args[1])
dontWriteHDF5(fds) <- TRUE

fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0, filter = FALSE)
fds <- fds[mcols(fds, type = "j")[, "passed"], ]

# FRASER 1.99.x/2.0 has optimHyperParams(); FRASER >= 2.4.2 has estimateBestQ()
# (useOHT=FALSE is the grid search of the older function).
ns <- asNamespace("FRASER")
if (exists("optimHyperParams", envir = ns, inherits = FALSE)) {
    search <- get("optimHyperParams", envir = ns)
} else {
    search <- function(fds, type, q_param, plot) {
        estimateBestQ(fds, type = type, useOHT = FALSE, q_param = q_param,
                      plot = plot)
    }
}

set.seed(42)
fds <- search(fds, type = "jaccard", q_param = q, plot = FALSE)

res <- tryCatch(hyperParams(fds, type = "jaccard", all = TRUE),
                error = function(e) metadata(fds)[["hyperParams_jaccard"]])
stopifnot(!is.null(res), nrow(res) >= 1)
write.table(as.data.frame(res), paste0(q, "_hyperparams.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
