### get project script
FD_REPO = "/hpc/group/igvf/kk319/repo"
FD_PRJ  = file.path(FD_REPO, "Proj_IGVF_BlueSTARR")
FD_EXE  = file.path(FD_PRJ,  "scripts")

### source
txt_fpath = file.path(FD_EXE, "config_project.R")
source(txt_fpath)