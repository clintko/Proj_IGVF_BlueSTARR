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
TXT_FPATH_MOTIF=$2
TXT_FPATH_OUT=$3
NUM_FLANK_LEFT=$4
NUM_BATCH_SIZE=$5

# =========================
# Execute
# -------------------------

### set script
TXT_FPATH_APP=./run_script.sh
TXT_FPATH_EXE=./run_motifdelta_01_scan.py

### run script
${TXT_FPATH_APP} python ${TXT_FPATH_EXE}  \
    --txt_fpath_fasta  "${TXT_FPATH_INP}" \
    --txt_fpath_motif  "${TXT_FPATH_MOTIF}" \
    --txt_fpath_output "${TXT_FPATH_OUT}"  \
    --num_flank_left   "${NUM_FLANK_LEFT}" \
    --batch_size       "${NUM_BATCH_SIZE}"
    
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

