#!/bin/bash
#SBATCH -J test_alphagenome_dnaclient
#SBATCH -A majoroslab                
#SBATCH -p igvf                      
#SBATCH --cpus-per-task=1 
#SBATCH --mem=2G                   
#SBATCH -o /hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/log/test.alphagenome.dnaclient.out 
#SBATCH -e /hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/log/test.alphagenome.dnaclient.err
#SBATCH --chdir=/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/scripts

set -euo pipefail
run_alphagenome python test_alphagenome_dnaclient.py 