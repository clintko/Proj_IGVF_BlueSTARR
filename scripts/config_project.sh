#!/bin/bash


#########################################
### Source files
### =====================================

### get path of the script
### Stackoverflow: how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

### source files
source ${SCRIPT_DIR}/config_path_duke_dcc_igvf.sh
source ${SCRIPT_DIR}/fun_utils.sh

########################################
### configuration for jupyter
### =====================================

### Hack to handle broken pipes - IGNORE.
### ex: suppress errors in `zcat ... | head`
cleanup () {
    :
}
trap "cleanup" SIGPIPE