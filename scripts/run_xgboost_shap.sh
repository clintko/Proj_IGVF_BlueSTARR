#!/bin/bash
set -euo pipefail

# =========================
# Environment
# -------------------------

timer_start=$(date +%s)
echo "Hostname:          " "$(hostname)"
echo "Slurm Array Index: " "${SLURM_ARRAY_TASK_ID-NA}"
echo "Time Stamp:        " "$(date +"%m-%d-%y+%T")"
echo "PWD:               " "$(pwd)"
echo

source fun_utils.sh

# =========================
# Defaults
# -------------------------
FP_INP_X=""
FP_INP_M=""
FP_OUT_PREFIX=""

usage() {
  cat <<EOF
Usage:
  $0 --fp_inp_x FP_INP_X --fp_inp_m FP_INP_M --fp_out_prefix FP_OUT_PREFIX

Required:
  --fp_inp_x       Path to X (.rds)
  --fp_inp_m       Path to model (.bin)
  --fp_out_prefix  Output prefix (full path; no extension)

Example:
  $0 \\
    --fp_inp_x /path/pilot.ori.X.rds \\
    --fp_inp_m /path/pilot.ori.xgb_model.bin \\
    --fp_out_prefix /path/out/variant_closed_gof_bluestarr.pilot.ori
EOF
}

# =========================
# Parse args
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fp_inp_x) FP_INP_X="$2"; shift 2 ;;
    --fp_inp_m) FP_INP_M="$2"; shift 2 ;;
    --fp_out_prefix) FP_OUT_PREFIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

# =========================
# Checks
# -------------------------
if [[ -z "${FP_INP_X}" ]];      then echo "ERROR: --fp_inp_x is required."; usage; exit 1; fi
if [[ -z "${FP_INP_M}" ]];      then echo "ERROR: --fp_inp_m is required."; usage; exit 1; fi
if [[ -z "${FP_OUT_PREFIX}" ]]; then echo "ERROR: --fp_out_prefix is required."; usage; exit 1; fi

if [[ ! -f "${FP_INP_X}" ]]; then echo "ERROR: FP_INP_X not found: ${FP_INP_X}"; exit 1; fi
if [[ ! -f "${FP_INP_M}" ]]; then echo "ERROR: FP_INP_M not found: ${FP_INP_M}"; exit 1; fi

FD_OUT="$(dirname "${FP_OUT_PREFIX}")"
mkdir -p "${FD_OUT}"

# =========================
# Show config
# -------------------------
echo "### FP_INP_X=${FP_INP_X}"
echo "### FP_INP_M=${FP_INP_M}"
echo "### FP_OUT_PREFIX=${FP_OUT_PREFIX}"
echo

# =========================
# Run
# -------------------------
TXT_FPATH_APP=./run_script.sh
TXT_FPATH_EXE=./run_xgboost_shap.R

if [[ ! -f "${TXT_FPATH_APP}" ]]; then echo "ERROR: app script not found: ${TXT_FPATH_APP}"; exit 1; fi
if [[ ! -f "${TXT_FPATH_EXE}" ]]; then echo "ERROR: R script not found: ${TXT_FPATH_EXE}"; exit 1; fi

${TXT_FPATH_APP} Rscript "${TXT_FPATH_EXE}" \
    --fp_inp_x "${FP_INP_X}" \
    --fp_inp_m "${FP_INP_M}" \
    --fp_out_prefix "${FP_OUT_PREFIX}"

# =========================
# Finish
# -------------------------
timer=$(date +%s)
runtime=$(( timer - timer_start ))
echo
echo "Done!"
echo "Run Time: $(displaytime ${runtime})"
echo