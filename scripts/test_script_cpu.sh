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

# ---- execute ----
echo "Env (Python)"
python  -c "import sys; print(sys.version)"
echo

echo "Env (R)"
Rscript -e 'print(R.version)'
echo

# ---- end message ----
timer=`date +%s`
runtime=$(( timer - timer_start ))
echo
echo 'Done!'
echo "Run Time: $(displaytime ${runtime})"
echo