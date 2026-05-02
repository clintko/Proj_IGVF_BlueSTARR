#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(Matrix)
  library(xgboost)
})

### ============================================================
### helpers
### ============================================================
fun_stop_if_not_file <- function(fp) {
  if (!file.exists(fp)) {
    stop(sprintf("ERROR: file not found: %s", fp))
  }
}

fun_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

### ============================================================
### parse args
### ============================================================
option_list <- list(
  make_option(c("--fp_inp_x"), type = "character", default = NA,
              help = "Path to X (.rds)."),
  make_option(c("--fp_inp_m"), type = "character", default = NA,
              help = "Path to xgboost model (.bin)."),
  make_option(c("--fp_out_prefix"), type = "character", default = NA,
              help = "Output prefix (full path; no extension).")
)

opt <- parse_args(OptionParser(option_list = option_list))

FP_INP_X      <- opt$fp_inp_x
FP_INP_M      <- opt$fp_inp_m
FP_OUT_PREFIX <- opt$fp_out_prefix

if (is.na(FP_INP_X)      || FP_INP_X      == "") stop("ERROR: --fp_inp_x is required.")
if (is.na(FP_INP_M)      || FP_INP_M      == "") stop("ERROR: --fp_inp_m is required.")
if (is.na(FP_OUT_PREFIX) || FP_OUT_PREFIX == "") stop("ERROR: --fp_out_prefix is required.")

fun_stop_if_not_file(FP_INP_X)
fun_stop_if_not_file(FP_INP_M)

### ============================================================
### outputs
### ============================================================
FP_OUT_SHAP_RDS      <- paste0(FP_OUT_PREFIX, ".shap.rds")
FP_OUT_SHAP_TSV      <- paste0(FP_OUT_PREFIX, ".shap.summary.tsv")
FP_OUT_SHAP_META_TXT <- paste0(FP_OUT_PREFIX, ".shap.meta.txt")

### ============================================================
### read inputs
### ============================================================
cat("### start time:", fun_now(), "\n")
cat("### FP_INP_X=", FP_INP_X, "\n", sep = "")
cat("### FP_INP_M=", FP_INP_M, "\n", sep = "")
cat("### FP_OUT_PREFIX=", FP_OUT_PREFIX, "\n", sep = "")
cat("\n")

cat("### read X\n")
mat_X <- readRDS(FP_INP_X)

### convert to sparse dgCMatrix if needed
if (inherits(mat_X, "dgCMatrix") || inherits(mat_X, "sparseMatrix")) {
  mat_X_sparse <- mat_X
} else if (is.matrix(mat_X)) {
  mat_X_sparse <- Matrix::Matrix(mat_X, sparse = TRUE)
} else if (is.data.frame(mat_X)) {
  mat_X_sparse <- Matrix::Matrix(as.matrix(mat_X), sparse = TRUE)
} else {
  stop(sprintf("Unsupported X type: %s", paste(class(mat_X), collapse = ", ")))
}

cat("### read model\n")
fit_xgb <- xgboost::xgb.load(FP_INP_M)

### build DMatrix (no label needed for SHAP)
cat("### build DMatrix\n")
dmat <- xgboost::xgb.DMatrix(data = mat_X_sparse)

### ============================================================
### compute SHAP
### ============================================================
cat("### compute SHAP via predict(predcontrib=TRUE)\n")
time_shap_start <- Sys.time()

### NOTE:
### - Returns a dense numeric matrix with columns = features + "BIAS"
### - SHAP values sum + BIAS = prediction for each row
mat_shap <- predict(fit_xgb, dmat, predcontrib = TRUE)

time_shap_end <- Sys.time()
dur_shap <- as.numeric(difftime(time_shap_end, time_shap_start, units = "mins"))

if (is.null(dim(mat_shap))) {
  stop("ERROR: SHAP output is not a matrix. Check xgboost predict() output.")
}

cat("### SHAP done\n")
cat("### shap_dim=", nrow(mat_shap), "x", ncol(mat_shap), "\n")
cat("### shap_minutes=", sprintf("%.4f", dur_shap), "\n")
cat("\n")

### ============================================================
### align feature names
### ============================================================
### xgboost returns "BIAS" as the last column (usually)
### try to set feature column names from X colnames (if present)
vec_feat <- colnames(mat_X_sparse)
if (!is.null(vec_feat) && length(vec_feat) + 1 == ncol(mat_shap)) {
  colnames(mat_shap) <- c(vec_feat, "BIAS")
}

### ============================================================
### save SHAP matrix
### ============================================================
cat("### save SHAP matrix (rds)\n")
saveRDS(mat_shap, FP_OUT_SHAP_RDS)

### ============================================================
### build SHAP summary table (exclude BIAS)
### ============================================================
cat("### build SHAP summary table\n")

if (ncol(mat_shap) < 2) {
  stop("ERROR: SHAP matrix has <2 columns; cannot separate BIAS and features.")
}

mat_shap_feat <- mat_shap[, 1:(ncol(mat_shap) - 1), drop = FALSE]

### per-feature summary
vec_mean_abs <- colMeans(abs(mat_shap_feat))
vec_mean     <- colMeans(mat_shap_feat)

### fraction nonzero (useful signal even though SHAP is dense-ish)
vec_frac_nz  <- colMeans(mat_shap_feat != 0)

dat_sum <- data.frame(
  feature = colnames(mat_shap_feat),
  mean_abs_shap = as.numeric(vec_mean_abs),
  mean_shap     = as.numeric(vec_mean),
  frac_nonzero  = as.numeric(vec_frac_nz),
  stringsAsFactors = FALSE
)

dat_sum <- dat_sum[order(dat_sum$mean_abs_shap, decreasing = TRUE), , drop = FALSE]
dat_sum$rank <- seq_len(nrow(dat_sum))

readr::write_tsv(dat_sum, FP_OUT_SHAP_TSV)

### ============================================================
### meta
### ============================================================
vec_meta <- c(
  paste0("start_time: ", fun_now()),
  paste0("FP_INP_X: ", FP_INP_X),
  paste0("FP_INP_M: ", FP_INP_M),
  paste0("FP_OUT_PREFIX: ", FP_OUT_PREFIX),
  paste0("nrow: ", nrow(mat_X_sparse)),
  paste0("ncol: ", ncol(mat_X_sparse)),
  paste0("shap_rows: ", nrow(mat_shap)),
  paste0("shap_cols: ", ncol(mat_shap)),
  paste0("shap_minutes: ", sprintf("%.4f", dur_shap)),
  paste0("xgboost_version: ", as.character(packageVersion("xgboost"))),
  paste0("Matrix_version: ", as.character(packageVersion("Matrix"))),
  paste0("R_version: ", R.version.string),
  paste0("note: shap matrix includes last column 'BIAS' if colnames were assigned.")
)

writeLines(vec_meta, FP_OUT_SHAP_META_TXT, useBytes = TRUE)

cat("### saved:\n")
cat("###   ", FP_OUT_SHAP_RDS, "\n", sep = "")
cat("###   ", FP_OUT_SHAP_TSV, "\n", sep = "")
cat("###   ", FP_OUT_SHAP_META_TXT, "\n", sep = "")
cat("\n")
cat("### done:", fun_now(), "\n")