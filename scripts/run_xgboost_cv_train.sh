#!/bin/bash
set -euo pipefail

# =========================
# Environment
# -------------------------

### start message
timer_start=`date +%s`
echo "Hostname:          " $(hostname)
echo "Slurm Array Index: " ${SLURM_ARRAY_TASK_ID-NA}
echo "Time Stamp:        " $(date +"%m-%d-%y+%T")
echo "PWD:               " $(pwd)
echo

### load helper function
source fun_utils.sh

### ----------------------------
### defaults
### ----------------------------
FP_INP_X=""
FP_INP_Y=""
FP_INP_PARAMS=""
FP_OUT_PREFIX=""

usage() {
  cat <<EOF
Usage:
  $0 --fp_inp_x FP_INP_X --fp_inp_y FP_INP_Y --fp_inp_params FP_INP_PARAMS --fp_out_prefix FP_OUT_PREFIX [options]

Required:
  --fp_inp_x           Path to X (.rds) (e.g., pilot.ori.X.rds or pilot.nuc.X.rds)
  --fp_inp_y           Path to y (.rds)
  --fp_inp_params      Path to params (.json)
  --fp_out_prefix  Output prefix (full path prefix; no extension)

Example:
  $0 \\
    --fp_inp_x /path/pilot.ori.X.rds \\
    --fp_inp_y /path/pilot.y.rds \\
    --fp_inp_params /path/params.json \\
    --fp_out_prefix /path/out/variant_closed_gof_bluestarr.pilot.ori

EOF
}

### ----------------------------
### parse args
### ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fp_inp_x) FP_INP_X="$2"; shift 2 ;;
    --fp_inp_y) FP_INP_Y="$2"; shift 2 ;;
    --fp_inp_params)  FP_INP_PARAMS="$2";  shift 2 ;;
    --fp_out_prefix)  FP_OUT_PREFIX="$2";  shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

### ----------------------------
### check args
### ----------------------------
if [[ -z "${FP_INP_X}" ]]; then
  echo "ERROR: --fp_inp_x is required."
  usage
  exit 1
fi

if [[ -z "${FP_INP_Y}" ]]; then
  echo "ERROR: --fp_inp_y is required."
  usage
  exit 1
fi

if [[ -z "${FP_INP_PARAMS}" ]]; then
  echo "ERROR: --fp_inp_params is required."
  usage
  exit 1
fi

if [[ -z "${FP_OUT_PREFIX}" ]]; then
  echo "ERROR: --fp_out_prefix is required."
  usage
  exit 1
fi

if [[ ! -f "${FP_INP_X}" ]]; then
  echo "ERROR: FP_INP_X not found: ${FP_INP_X}"
  exit 1
fi

if [[ ! -f "${FP_INP_Y}" ]]; then
  echo "ERROR: FP_INP_Y not found: ${FP_INP_Y}"
  exit 1
fi

if [[ ! -f "${FP_INP_PARAMS}" ]]; then
  echo "ERROR: FP_INP_PARAMS not found: ${FP_INP_PARAMS}"
  exit 1
fi

### ensure output directory exists
FD_OUT="$(dirname "${FP_OUT_PREFIX}")"
mkdir -p "${FD_OUT}"

### ----------------------------
### show config
### ----------------------------
echo "### FP_INP_X=${FP_INP_X}"
echo "### FP_INP_Y=${FP_INP_Y}"
echo "### FP_INP_PARAMS=${FP_INP_PARAMS}"
echo "### FP_OUT_PREFIX=${FP_OUT_PREFIX}"
echo

### ----------------------------
### run
### ----------------------------

### set script
TXT_FPATH_APP=./run_script.sh
TXT_FPATH_EXE=./run_xgboost_cv_train.R

### double check the script exist
if [[ ! -f "${TXT_FPATH_APP}" ]]; then
  echo "ERROR: app script not found: ${TXT_FPATH_APP}"
  exit 1
fi

if [[ ! -f "${TXT_FPATH_EXE}" ]]; then
  echo "ERROR: R script not found: ${TXT_FPATH_EXE}"
  exit 1
fi

### run script
${TXT_FPATH_APP} Rscript "${TXT_FPATH_EXE}" \
  --fp_inp_x "${FP_INP_X}" \
  --fp_inp_y "${FP_INP_Y}" \
  --fp_inp_params "${FP_INP_PARAMS}" \
  --fp_out_prefix "${FP_OUT_PREFIX}"

# =========================
# Finish
# -------------------------

### print end message
timer=$(date +%s)
runtime=$(( timer - timer_start ))
echo
echo 'Done!'
echo "Run Time: $(displaytime ${runtime})"
echo