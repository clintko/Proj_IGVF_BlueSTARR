#!/bin/bash

### Helper function to display script running time
### Reference:
###     https://unix.stackexchange.com/questions/27013/displaying-seconds-as-days-hours-mins-seconds
function displaytime {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d days ' $D
  (( $H > 0 )) && printf '%d hours ' $H
  (( $M > 0 )) && printf '%d minutes ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d seconds\n' $S
}

### Helper function: wrapper of cat for handling multiple file types
function fun_cat {
    FPATH=$1
    
    ### check if file is a symlink, if yes, get the file path
    if [[ -h ${FPATH} ]]; then
        FPATH=$(readlink -f ${FPATH})
    fi
    
    ### check if file is compress, if yes: `zcat`, else: `cat`
    if (file ${FPATH} | grep -q compressed); then
        zcat ${FPATH}
    else
        cat  ${FPATH}
    fi
}