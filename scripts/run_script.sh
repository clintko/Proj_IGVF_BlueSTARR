#!/bin/bash

#########################################
### Wrapper of project container
### =====================================

### container image
DIR="/hpc/group/igvf/kk319/container/project"
IMG="singularity_proj_igvf_bluestarr.sif"
APP="${DIR}/${IMG}"

### point PYTHONPATH to my scripts & package directory
FD_EXE="/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/scripts"
FD_MOTIFDELTA="/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/scripts/motifdelta/src"
export PYTHONPATH="${FD_MOTIFDELTA}:${FD_EXE}:${PYTHONPATH:-}"

### get command and arguments
### stored in an array to preserve each argument exactly as-is
CMD=("$@")

### execute the command
singularity exec \
    --env NUMBA_CACHE_DIR=/tmp \
    --env PYTHONPATH="${PYTHONPATH}" \
    -B "${PWD}:${PWD}" \
    -B /tmp:/tmp \
    -B /hpc:/hpc \
    -B /datacommons:/datacommons \
    "${APP}" "${CMD[@]}"
