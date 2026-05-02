#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
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

fun_write_lines <- function(x, fp) {
  writeLines(text = x, con = fp, useBytes = TRUE)
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
  make_option(c("--fp_inp_y"), type = "character", default = NA,
              help = "Path to y (.rds)."),
  make_option(c("--fp_inp_params"), type = "character", default = NA,
              help = "Path to params (.json)."),
  make_option(c("--fp_out_prefix"), type = "character", default = NA,
              help = "Output prefix (full path; no extension).")
)

### parse parameters
opt <- parse_args(OptionParser(option_list = option_list))

FP_INP_X      <- opt$fp_inp_x
FP_INP_Y      <- opt$fp_inp_y
FP_INP_PARAMS <- opt$fp_inp_params
FP_OUT_PREFIX <- opt$fp_out_prefix

### ensure parameters not null or NA
if (is.na(FP_INP_X)      || FP_INP_X      == "") stop("ERROR: --fp_inp_x is required.")
if (is.na(FP_INP_Y)      || FP_INP_Y      == "") stop("ERROR: --fp_inp_y is required.")
if (is.na(FP_INP_PARAMS) || FP_INP_PARAMS == "") stop("ERROR: --fp_inp_params is required.")
if (is.na(FP_OUT_PREFIX) || FP_OUT_PREFIX == "") stop("ERROR: --fp_out_prefix is required.")

### ensure input file exist
fun_stop_if_not_file(FP_INP_X)
fun_stop_if_not_file(FP_INP_Y)
fun_stop_if_not_file(FP_INP_PARAMS)

### ============================================================
### outputs
### ============================================================
FP_OUT_CV_RDS        <- paste0(FP_OUT_PREFIX, ".cv.rds")
FP_OUT_EVAL_LOG_TSVG <- paste0(FP_OUT_PREFIX, ".cv.evaluation_log.tsv")
FP_OUT_MODEL_BIN     <- paste0(FP_OUT_PREFIX, ".xgb_model.bin")
FP_OUT_SUMMARY_TXT   <- paste0(FP_OUT_PREFIX, ".summary.txt")
FP_OUT_PARAMS_USED   <- paste0(FP_OUT_PREFIX, ".params.used.json")

FP_OUT_YHAT_RDS      <- paste0(FP_OUT_PREFIX, ".yhat.rds")
FP_OUT_YHAT_TSV      <- paste0(FP_OUT_PREFIX, ".yhat.tsv")
FP_OUT_RESID_RDS     <- paste0(FP_OUT_PREFIX, ".resid.rds")
FP_OUT_METRICS_TXT   <- paste0(FP_OUT_PREFIX, ".train_metrics.txt")

### ============================================================
### read params json
### ============================================================
cfg <- jsonlite::fromJSON(FP_INP_PARAMS)

### required keys
if (is.null(cfg$seed))  stop("ERROR: params json missing: seed")
if (is.null(cfg$cv))    stop("ERROR: params json missing: cv")
if (is.null(cfg$model)) stop("ERROR: params json missing: model")

if (is.null(cfg$cv$nfold))                 stop("ERROR: params json missing: cv$nfold")
if (is.null(cfg$cv$nrounds))               stop("ERROR: params json missing: cv$nrounds")
if (is.null(cfg$cv$early_stopping_rounds)) stop("ERROR: params json missing: cv$early_stopping_rounds")

### optional keys
if (is.null(cfg$cv$verbose)) cfg$cv$verbose <- 0L
if (is.null(cfg$train) || is.null(cfg$train$verbose)) {
    cfg$train <- list(verbose = 1L)
}

### normalize types
NUM_SEED          <- as.integer(cfg$seed)
NUM_NFOLD         <- as.integer(cfg$cv$nfold)
NUM_NROUNDS       <- as.integer(cfg$cv$nrounds)
NUM_EARLY_STOP    <- as.integer(cfg$cv$early_stopping_rounds)
NUM_VERBOSE_CV    <- as.integer(cfg$cv$verbose)
NUM_VERBOSE_TRAIN <- as.integer(cfg$train$verbose)

### model parameters passed to xgboost
lst_param <- cfg$model

### print summary (useful for SLURM logs)
cat("### parameters loaded\n")
cat("NUM_SEED:", NUM_SEED, "\n")
cat("NUM_NFOLD:", NUM_NFOLD, "\n")
cat("NUM_NROUNDS:", NUM_NROUNDS, "\n")
cat("NUM_EARLY_STOP:", NUM_EARLY_STOP, "\n")
cat("NUM_VERBOSE_CV:", NUM_VERBOSE_CV, "\n")
cat("NUM_VERBOSE_TRAIN:", NUM_VERBOSE_TRAIN, "\n")
cat("model parameters:\n")
print(lst_param)
cat("\n")

### ensure lst_param does not contain verbose
if (!is.null(lst_param$verbose)) lst_param$verbose <- NULL

### ============================================================
### read data
### ============================================================
cat("### start time:",    fun_now(),     "\n")
cat("### FP_INP_X=",      FP_INP_X,      "\n", sep = "")
cat("### FP_INP_Y=",      FP_INP_Y,      "\n", sep = "")
cat("### FP_INP_PARAMS=", FP_INP_PARAMS, "\n", sep = "")
cat("### FP_OUT_PREFIX=", FP_OUT_PREFIX, "\n", sep = "")
cat("\n")

cat("### read X\n")
mat_X <- readRDS(FP_INP_X)

cat("### read y\n")
vec_y <- readRDS(FP_INP_Y)

### basic checks
if (is.null(dim(mat_X))) {
    stop("ERROR: X must be a matrix/data.frame-like object with dim().")
}
if (length(vec_y) != nrow(mat_X)) {
    stop(sprintf("ERROR: length(y) != nrow(X): %d vs %d", length(vec_y), nrow(mat_X)))
}

### ============================================================
### make sparse + DMatrix
### ============================================================
cat("### convert X to sparse dgCMatrix\n")

### convert X to sparse dgCMatrix safely
if (inherits(mat_X, "dgCMatrix") || inherits(mat_X, "sparseMatrix")) {
    ### already sparse
    mat_X_sparse <- mat_X

} else if (is.matrix(mat_X)) {
    ### dense matrix → sparse
    mat_X_sparse <- Matrix::Matrix(mat_X, sparse = TRUE)

} else if (is.data.frame(mat_X)) {
    ### data.frame → matrix → sparse
    mat_X_sparse <- Matrix::Matrix(as.matrix(mat_X), sparse = TRUE)

} else {
    stop(sprintf("Unsupported X type: %s", paste(class(mat_X), collapse = ", ")))
}

### free dense object if possible
rm(mat_X)
gc()

### combine X & y for later call in xgboost modeling
cat("### build DMatrix\n")
dtrain <- xgboost::xgb.DMatrix(
  data  = mat_X_sparse,
  label = vec_y
)

### ============================================================
### run CV
### ============================================================
set.seed(NUM_SEED)

cat("### xgb.cv\n")
time_cv_start <- Sys.time()

fit_xgb_cv <- xgboost::xgb.cv(
  params                = lst_param,
  data                  = dtrain,
  nrounds               = NUM_NROUNDS,
  nfold                 = NUM_NFOLD,
  early_stopping_rounds = NUM_EARLY_STOP,
  verbose               = NUM_VERBOSE_CV
)

time_cv_end <- Sys.time()
dur_cv <- as.numeric(difftime(time_cv_end, time_cv_start, units = "mins"))

### ------------------------------------------------------------
### choose best iteration from CV
### ------------------------------------------------------------

### extract evaluation log
eval_log <- fit_xgb_cv$evaluation_log

### extract best iteration reported by xgboost
best_iter <- fit_xgb_cv$best_iteration
cat("### best_iter (raw): ",
    ifelse(is.null(best_iter), "NULL", as.character(best_iter)),
    "\n", sep = "")

### fallback logic
### note: in some xgboost builds, best_iteration may not be populated
### even when early stopping occurred
if (is.null(best_iter) || is.na(best_iter)) {

    ### ensure evaluation log exists
    if (is.null(eval_log) || nrow(eval_log) == 0) {
        stop("ERROR: evaluation_log is empty; cannot determine best_iter.")
    }

    ### check if CV ended before NUM_NROUNDS
    ### this typically indicates early stopping triggered
    if (nrow(eval_log) < NUM_NROUNDS) {
        cat("### CV ended early (eval_log rows < NUM_NROUNDS); likely early stopping.\n")
    }

    ### determine best iteration from evaluation log
    ### (minimum test_rmse_mean)
    best_iter <- which.min(eval_log$test_rmse_mean)

    cat("### WARNING: best_iteration is NULL/NA.\n")
    cat("### fallback best_iter (argmin test_rmse_mean): ",
        best_iter, "\n", sep = "")
}

### extract best rmse
best_test_rmse <- eval_log$test_rmse_mean[best_iter]

### compute early stopping rounds
### (number of rounds after best_iter before CV stopped)
num_early_stop_round <- nrow(eval_log) - best_iter

### report CV results
cat("\n")
cat("### CV done\n")
cat("### best_iter=",        best_iter, "\n", sep = "")
cat("### best_test_rmse=",   sprintf("%.8f", best_test_rmse), "\n", sep = "")
cat("### eval_log rows: ",   nrow(eval_log), "\n", sep = "")
cat("### early_stop_round=", num_early_stop_round, "\n", sep = "")
cat("### cv_minutes=",       sprintf("%.4f", dur_cv), "\n", sep = "")
cat("\n")

### ============================================================
### train final model
### ============================================================
set.seed(NUM_SEED)

cat("### xgb.train (final)\n")
time_train_start <- Sys.time()

fit_xgb_final <- xgboost::xgb.train(
  params  = lst_param,
  data    = dtrain,
  nrounds = best_iter,
  verbose = NUM_VERBOSE_TRAIN
)

time_train_end <- Sys.time()
dur_train <- as.numeric(difftime(time_train_end, time_train_start, units = "mins"))

cat("\n")
cat("### train done\n")
cat("### train_minutes=", sprintf("%.4f", dur_train), "\n", sep = "")
cat("\n")

### ============================================================
### predict on training X
### ============================================================
cat("### predict (train)\n")
time_pred_start <- Sys.time()

vec_yhat <- predict(fit_xgb_final, dtrain)

time_pred_end <- Sys.time()
dur_pred <- as.numeric(difftime(time_pred_end, time_pred_start, units = "mins"))

### basic sanity check
if (length(vec_yhat) != length(vec_y)) {
    stop(sprintf("ERROR: length(yhat) != length(y): %d vs %d",
                 length(vec_yhat), length(vec_y)))
}

### residuals + simple metrics
vec_resid <- vec_yhat - vec_y
mse  <- mean((vec_resid)^2)
rmse <- sqrt(mse)
mae  <- mean(abs(vec_resid))
r2   <- 1 - sum((vec_resid)^2) / sum((vec_y - mean(vec_y))^2)

cat("### predict done\n")
cat("### pred_minutes=", sprintf("%.4f", dur_pred), "\n", sep = "")
cat("### train_rmse=",   sprintf("%.8f", rmse), "\n", sep = "")
cat("### train_mae=",    sprintf("%.8f", mae), "\n", sep = "")
cat("### train_r2=",     sprintf("%.8f", r2), "\n", sep = "")
cat("\n")


### ============================================================
### save outputs
### ============================================================

### save modeling results
cat("### save outputs: CV training and evaluation log\n")
saveRDS(fit_xgb_cv, FP_OUT_CV_RDS)
readr::write_tsv(eval_log, FP_OUT_EVAL_LOG_TSVG)

cat("### save outputs: Final model\n")
xgboost::xgb.save(fit_xgb_final, FP_OUT_MODEL_BIN)

cat("### save outputs: Model prediction\n")
saveRDS(vec_yhat,  FP_OUT_YHAT_RDS)
saveRDS(vec_resid, FP_OUT_RESID_RDS)

### ensure row ids exist
vec_id <- rownames(mat_X_sparse)
if (is.null(vec_id)) {
    vec_id <- as.character(seq_len(nrow(mat_X_sparse)))
}

dat_yhat <- data.frame(
    id = vec_id,
    y  = as.numeric(vec_y),
    yhat  = as.numeric(vec_yhat),
    resid = as.numeric(vec_resid),
    stringsAsFactors = FALSE
)
readr::write_tsv(dat_yhat, FP_OUT_YHAT_TSV)

### small metrics text
vec_metrics <- c(
    paste0("train_rmse: ", sprintf("%.8f", rmse)),
    paste0("train_mae: ",  sprintf("%.8f", mae)),
    paste0("train_r2: ",   sprintf("%.8f", r2)),
    paste0("pred_minutes: ", sprintf("%.4f", dur_pred))
)
fun_write_lines(vec_metrics, FP_OUT_METRICS_TXT)

### save exact params used (including derived numbers)
cat("### save outputs: Exact params used\n")

lst_params_used <- list(
    fp_inp_x = FP_INP_X,
    fp_inp_y = FP_INP_Y,
    fp_inp_params = FP_INP_PARAMS,
    fp_out_prefix = FP_OUT_PREFIX,
    seed = NUM_SEED,
    cv = list(
        nfold = NUM_NFOLD,
        nrounds = NUM_NROUNDS,
        early_stopping_rounds = NUM_EARLY_STOP,
        verbose = NUM_VERBOSE_CV
    ),
    train = list(
        verbose = NUM_VERBOSE_TRAIN,
        nrounds = best_iter
    ),
    model_params = lst_param,
    results = list(
        best_iter = best_iter,
        best_test_rmse = best_test_rmse,
        cv_minutes = dur_cv,
        train_minutes = dur_train,
        pred_minutes = dur_pred,
        train_rmse = rmse,
        train_mae = mae,
        train_r2 = r2
    )
)

jsonlite::write_json(
  lst_params_used,
  FP_OUT_PARAMS_USED,
  auto_unbox = TRUE,
  pretty = TRUE
)

### summary text (human readable)
cat("### save outputs: summary text of training\n")

vec_summary <- c(
  paste0("start_time: ", fun_now()),
  paste0("FP_INP_X: ", FP_INP_X),
  paste0("FP_INP_Y: ", FP_INP_Y),
  paste0("FP_INP_PARAMS: ", FP_INP_PARAMS),
  paste0("FP_OUT_PREFIX: ", FP_OUT_PREFIX),
  "",
  paste0("nrow: ", nrow(mat_X_sparse)),
  paste0("ncol: ", ncol(mat_X_sparse)),
  "",
  paste0("best_iter: ", best_iter),
  paste0("best_test_rmse: ", sprintf("%.8f", best_test_rmse)),
  paste0("cv_minutes: ", sprintf("%.4f", dur_cv)),
  paste0("train_minutes: ", sprintf("%.4f", dur_train)),
  "",
  paste0("xgboost_version: ", as.character(packageVersion("xgboost"))),
  paste0("Matrix_version: ", as.character(packageVersion("Matrix"))),
  paste0("R_version: ", R.version.string)
)

fun_write_lines(vec_summary, FP_OUT_SUMMARY_TXT)

cat("### saved:\n")
cat("###   ", FP_OUT_CV_RDS, "\n", sep = "")
cat("###   ", FP_OUT_EVAL_LOG_TSVG, "\n", sep = "")
cat("###   ", FP_OUT_MODEL_BIN, "\n", sep = "")
cat("###   ", FP_OUT_PARAMS_USED, "\n", sep = "")
cat("###   ", FP_OUT_SUMMARY_TXT, "\n", sep = "")
cat("###   ", FP_OUT_YHAT_RDS, "\n", sep = "")
cat("###   ", FP_OUT_YHAT_TSV, "\n", sep = "")
cat("###   ", FP_OUT_RESID_RDS, "\n", sep = "")
cat("###   ", FP_OUT_METRICS_TXT, "\n", sep = "")
cat("\n")
cat("### done:", fun_now(), "\n")