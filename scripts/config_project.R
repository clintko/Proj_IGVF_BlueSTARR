# =====================================
# set project information
# =====================================

PROJECT <- "IGVF BlueSTARR"

# =====================================
# import common libraries
# =====================================

### import library: basic
library("tidyverse")

### import library: visualization
library("RColorBrewer")
library("cowplot")

# =====================================
# set other folders
# =====================================
FD_MLAB       <- "/hpc/group/majoroslab"
FD_IGVF       <- "/hpc/group/igvf"
FD_IGVF_KUEI  <- "/hpc/group/igvf/kk319"
FD_IGVF_REVA  <- "/hpc/group/igvf/revathy"

# =====================================
# set base working folder
# =====================================
FD_BASE <- "/hpc/group/igvf/kk319"
FD_REPO <- file.path(FD_BASE, "repo")
FD_WORK <- file.path(FD_BASE, "work")
FD_DATA <- file.path(FD_BASE, "data")
FD_ENVS <- file.path(FD_BASE, "envs")
FD_BINS <- file.path(FD_BASE, "bins")
FD_TEMP <- file.path(FD_BASE, "temp")
FD_SING <- file.path(FD_BASE, "container")

# =====================================
# set project folder
# =====================================
FD_PRJ <- file.path(FD_REPO, "Proj_IGVF_BlueSTARR")

FD_RES <- file.path(FD_PRJ, "results")
FD_EXE <- file.path(FD_PRJ, "scripts")
FD_DAT <- file.path(FD_PRJ, "data")
FD_NBK <- file.path(FD_PRJ, "notebooks")
FD_DOC <- file.path(FD_PRJ, "docs")
FD_LOG <- file.path(FD_PRJ, "log")
FD_REF <- file.path(FD_PRJ, "references")

# =====================================
# helper function
# =====================================
show_env <- function() {
  cat("BASE DIRECTORY (FD_BASE):", FD_BASE, "\n")
  cat("REPO DIRECTORY (FD_REPO):", FD_REPO, "\n")
  cat("WORK DIRECTORY (FD_WORK):", FD_WORK, "\n")
  cat("DATA DIRECTORY (FD_DATA):", FD_DATA, "\n\n")
  
  cat("You are working with     ", PROJECT, "\n")
  cat("PATH OF PROJECT (FD_PRJ):", FD_PRJ, "\n")
  cat("PROJECT RESULTS (FD_RES):", FD_RES, "\n")
  cat("PROJECT SCRIPTS (FD_EXE):", FD_EXE, "\n")
  cat("PROJECT DATA    (FD_DAT):", FD_DAT, "\n")
  cat("PROJECT NOTE    (FD_NBK):", FD_NBK, "\n")
  cat("PROJECT DOCS    (FD_DOC):", FD_DOC, "\n")
  cat("PROJECT LOG     (FD_LOG):", FD_LOG, "\n")
  cat("PROJECT REF     (FD_REF):", FD_REF, "\n")
  cat("\n")
}

### helper function to show table
fun_display_table = function(dat){
    dat %>%
        kableExtra::kable("html") %>%
        as.character() %>%
        IRdisplay::display_html()
}

### helper function to show table
fun_markdown_table = function(dat){
    dat %>% kableExtra::kable("markdown")    
}

