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

TXT_FPATH_INP=$1
TXT_FPATH_OUT=$2
NUM_POS0=$3
NUM_SEED=$4

# =========================
# Execute
# -------------------------

### set script
TXT_FPATH_APP=./run_script.sh
TXT_FPATH_EXE="./run_sequence_dinuc_shuffle.py"

### run script
${TXT_FPATH_APP} python ${TXT_FPATH_EXE} \
    --txt_finp "${TXT_FPATH_INP}" \
    --txt_fout "${TXT_FPATH_OUT}" \
    --num_pos0 "${NUM_POS0}" \
    --num_seed "${NUM_SEED}"

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

