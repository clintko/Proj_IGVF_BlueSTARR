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

TXT_FPATH_DELTA=$1
TXT_FPATH_PREFIX=$2

# =========================
# Execute
# -------------------------

### set script
TXT_FPATH_APP=./run_script.sh
TXT_FPATH_EXE=./run_motifdelta_03_merge.py

### run script
${TXT_FPATH_APP} python ${TXT_FPATH_EXE} \
    --txt_fpath_dir "${TXT_FPATH_DELTA}" \
    --txt_fpath_output_prefix "${TXT_FPATH_PREFIX}"

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

