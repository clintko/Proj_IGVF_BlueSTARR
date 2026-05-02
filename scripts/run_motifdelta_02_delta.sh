#!/bin/bash

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

# =========================
# I/O
# -------------------------

TXT_FPATH_SCAN=$1
TXT_FPATH_TBIND=$2
TXT_FPATH_OUT_PREFIX=$3

# =========================
# Execute
# -------------------------

### set script
TXT_FPATH_APP=./run_script.sh
TXT_FPATH_EXE=./run_motifdelta_02_delta.py

### run script
${TXT_FPATH_APP} python ${TXT_FPATH_EXE} \
    --txt_fpath_scan "${TXT_FPATH_SCAN}" \
    --txt_fpath_tbind "${TXT_FPATH_TBIND}" \
    --txt_fpath_output_prefix "${TXT_FPATH_OUT_PREFIX}"

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

