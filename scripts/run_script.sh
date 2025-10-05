#!/bin/bash

#########################################
### Wrapper of project container
### =====================================

### container image
DIR="/hpc/group/igvf/kk319/container/project"
IMG="singularity_proj_igvf_bluestarr.sif"
APP="${DIR}/${IMG}"

# point PYTHONPATH to my scripts/ directory
FD_EXE="/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/scripts"
export PYTHONPATH="${FD_EXE}:${PYTHONPATH}"

### get command and arguments
### stored in an array to preserve each argument exactly as-is
CMD=("$@")

### execute the command
singularity exec \
    --env NUMBA_CACHE_DIR=/tmp \
    --env PYTHONPATH="${PYTHONPATH}" \
    -B ${PWD} \
    -B /hpc:/hpc \
    -B /datacommons:/datacommons \
    ${APP} "${CMD[@]}"
