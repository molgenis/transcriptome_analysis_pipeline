#' DESeq2 differential expression between the groups of the samplesheet
#' Argument 1= merged gene counts (GeneID + one column per sample)
#' Argument 2= samplesheet (needs sampleID and group)
#' Argument 3= gene annotation table (gene_annotation.tsv, made from the GTF)
#' Argument 4= reference group (denominator); empty = all pairwise comparisons
#' Argument 5= covariates: comma separated samplesheet columns to correct for
#' Argument 6= minimum total count of a gene over all samples
#' Argument 7= number of CPUs
#'
#' Design: ~ covariates + group. For every comparison B_vs_A a positive
#' log2FoldChange means higher expression in B than in A.
#' Output: deseq2_<B>_vs_<A>.tsv, deseq2_summary.tsv, deseq2_normalized_counts.tsv,
#' deseq2_sample_info.tsv and deseq2_pca.pdf

suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
sheet_file  <- args[2]
anno_file   <- args[3]
reference   <- trimws(args[4])
covariates  <- trimws(strsplit(args[5], ",")[[1]])
covariates  <- covariates[covariates != ""]
min_count   <- as.numeric(args[6])
ncpu        <- as.integer(args[7])

bp <- if (ncpu > 1 && .Platform$OS.type == "unix") MulticoreParam(ncpu) else SerialParam()

# ---- counts and samples ---------------------------------------------------
cts <- fread(counts_file)
mat <- as.matrix(cts[, -1, with = FALSE])
rownames(mat) <- cts[[1]]

sheet <- fread(sheet_file)
if (!"group" %in% names(sheet)) stop("The samplesheet has no 'group' column")
sheet[, group := as.character(group)]
sheet <- sheet[!is.na(group) & !group %in% c("", "NA")]

missing <- setdiff(sheet$sampleID, colnames(mat))
if (length(missing) > 0) {
    message("Samples without counts are left out: ", paste(missing, collapse = ", "))
    sheet <- sheet[sampleID %in% colnames(mat)]
}

# groups need replicates
n_per_group <- table(sheet$group)
small <- names(n_per_group)[n_per_group < 2]
if (length(small) > 0) {
    message("Groups with a single sample are left out: ", paste(small, collapse = ", "))
    sheet <- sheet[!group %in% small]
}
if (length(unique(sheet$group)) < 2) {
    stop("DESeq2 needs at least two groups with at least two samples each")
}

lev <- sort(unique(sheet$group))
if (nzchar(reference)) {
    if (!reference %in% lev) {
        stop("Reference group '", reference, "' is not among the usable groups: ",
             paste(lev, collapse = ", "))
    }
    lev <- c(reference, setdiff(lev, reference))
}

# ---- covariates -----------------------------------------------------------
col_data <- data.frame(row.names = sheet$sampleID,
                       group = factor(sheet$group, levels = lev))
for (cv in covariates) {
    if (!cv %in% names(sheet)) stop("Covariate '", cv, "' is not a samplesheet column")
    v <- sheet[[cv]]
    if (anyNA(v) || any(as.character(v) %in% c("", "NA"))) {
        message("Covariate '", cv, "' has missing values and is not used")
        next
    }
    if (length(unique(v)) < 2) {
        message("Covariate '", cv, "' has a single value and is not used")
        next
    }
    col_data[[cv]] <- if (is.numeric(v)) v else factor(v)
}
used_cov <- setdiff(names(col_data), "group")
design <- as.formula(paste("~", paste(c(used_cov, "group"), collapse = " + ")))
message("Design: ", deparse(design))

# ---- DESeq2 ---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(mat[, rownames(col_data), drop = FALSE]),
                              colData = col_data, design = design)
dds <- dds[rowSums(counts(dds)) >= min_count, ]
message(nrow(dds), " genes with at least ", min_count, " counts in total")
dds <- DESeq(dds, parallel = ncpu > 1, BPPARAM = bp)

# ---- gene annotation (from the GTF) -----------------------------------------
anno <- fread(anno_file)
ids <- rownames(dds)
idx <- match(ids, anno$gene_id)
nover <- is.na(idx)
idx[nover] <- match(sub("\\..*$", "", ids[nover]), anno$gene_id_nover)
symbol <- anno$gene_name[idx]
symbol[!is.na(symbol) & symbol == ""] <- NA

# ---- contrasts --------------------------------------------------------------
lev <- levels(col_data$group)
pairs <- if (nzchar(reference)) {
    lapply(lev[-1], function(b) c(b, lev[1]))
} else {
    unlist(lapply(seq_len(length(lev) - 1), function(i)
        lapply((i + 1):length(lev), function(j) c(lev[j], lev[i]))), recursive = FALSE)
}

safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
have_ashr <- requireNamespace("ashr", quietly = TRUE)
summary_rows <- list()

for (pr in pairs) {
    b <- pr[1]; a <- pr[2]
    res <- results(dds, contrast = c("group", b, a), alpha = 0.05, parallel = ncpu > 1,
                   BPPARAM = bp)
    out <- data.table(EnsemblID = rownames(res), hgncSymbol = symbol,
                      baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
                      lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue,
                      padj = res$padj)
    if (have_ashr) {
        shr <- tryCatch(lfcShrink(dds, contrast = c("group", b, a), res = res,
                                  type = "ashr", quiet = TRUE),
                        error = function(e) NULL)
        if (!is.null(shr)) out[, log2FoldChangeShrunk := shr$log2FoldChange]
    }
    out <- out[order(padj, pvalue)]
    fwrite(out, sprintf("deseq2_%s_vs_%s.tsv", safe(b), safe(a)), sep = "\t")
    summary_rows[[length(summary_rows) + 1]] <- data.table(
        comparison = sprintf("%s_vs_%s", b, a),
        n_group_B = sum(col_data$group == b), n_group_A = sum(col_data$group == a),
        genes_tested = sum(!is.na(res$pvalue)),
        padj_0.05 = sum(res$padj < 0.05, na.rm = TRUE),
        up_in_B = sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
        down_in_B = sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE))
}
fwrite(rbindlist(summary_rows), "deseq2_summary.tsv", sep = "\t")

# ---- normalized counts, sample info, PCA ------------------------------------
norm <- as.data.table(counts(dds, normalized = TRUE), keep.rownames = "EnsemblID")
norm[, hgncSymbol := symbol]
setcolorder(norm, c("EnsemblID", "hgncSymbol"))
fwrite(norm, "deseq2_normalized_counts.tsv", sep = "\t")

info <- data.table(sampleID = colnames(dds), sizeFactor = sizeFactors(dds),
                   as.data.table(as.data.frame(colData(dds))[, setdiff(names(colData(dds)), "sizeFactor"), drop = FALSE]))
fwrite(info, "deseq2_sample_info.tsv", sep = "\t")

tryCatch({
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    m <- assay(vsd)
    top <- head(order(apply(m, 1, var), decreasing = TRUE), 500)
    pc <- prcomp(t(m[top, ]))
    pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
    cols <- as.integer(col_data$group)
    pdf("deseq2_pca.pdf", width = 7, height = 6)
    plot(pc$x[, 1], pc$x[, 2], col = cols, pch = 19,
         xlab = sprintf("PC1 (%s%%)", pct[1]), ylab = sprintf("PC2 (%s%%)", pct[2]),
         main = "PCA of the 500 most variable genes (VST)")
    legend("topright", legend = levels(col_data$group),
           col = seq_along(levels(col_data$group)), pch = 19, bty = "n")
    dev.off()
}, error = function(e) message("PCA plot skipped: ", conditionMessage(e)))
