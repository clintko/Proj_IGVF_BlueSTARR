#!/bin/bash

# ---- start message ----
timer_start=`date +%s`
echo "Hostname:          " $(hostname)
echo "Slurm Array Index: " ${SLURM_ARRAY_TASK_ID-NA}
echo "Time Stamp:        " $(date +"%m-%d-%y+%T")
echo

# ---- setup env ----
# load helper functions
source fun_utils.sh

# parse argument
# Example: 
#   run_alphagenome_variant_prediction.sh ${INPUT_FILE} ${OUTPUT_FILE} "2KB" "RNA_SEQ" "K562"
# For more information, run 
#   python run_alphagenome_variant_prediction.py -h
FP_INP=${1}
FP_OUT=${2}
TXT_SEQ_LEN=${3}
TXT_SCORER=${4}
TXT_BIOSAMPLE=${5}

# ---- show input ----
echo "Input: " ${FP_INP}
echo
echo "show first few lines of input"
fun_cat ${FP_INP} | head -n 3
echo

# ---- execute ----
run_alphagenome python run_alphagenome_variant_prediction.py \
    -i ${FP_INP} \
    -o ${FP_OUT} \
    -s ${TXT_SEQ_LEN} \
    -c ${TXT_SCORER} \
    -b ${TXT_BIOSAMPLE} \
    --min-delay 0.3
    
# ---- show input ----
echo
echo "Output: " ${FP_OUT}
echo
echo "show first few lines of output:"
fun_cat ${FP_OUT} | head -n 3
echo

# ---- end message ----
timer=`date +%s`
runtime=$(( timer - timer_start ))
echo
echo 'Done!'
echo "Run Time: $(displaytime ${runtime})"
echo