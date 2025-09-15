import sys, os

PROJECT = "IGVF BlueSTARR"

# =====================================
# set other folders
# =====================================
FD_MLAB = "/hpc/group/majoroslab"
FD_IGVF = "/hpc/group/igvf"
FD_IGVF_KUEI = "/hpc/group/igvf/kk319"
FD_IGVF_REVA = "/hpc/group/igvf/revathy"

# =====================================
# set base working folder
# =====================================
FD_BASE = "/hpc/group/igvf/kk319"
FD_REPO = f"{FD_BASE}/repo"
FD_WORK = f"{FD_BASE}/work"
FD_DATA = f"{FD_BASE}/data"
FD_ENVS = f"{FD_BASE}/envs"
FD_BINS = f"{FD_BASE}/bins"
FD_TEMP = f"{FD_BASE}/temp"
FD_SING = f"{FD_BASE}/container"

# =====================================
# set project folder
# =====================================
FD_PRJ = f"{FD_REPO}/Proj_IGVF_BlueSTARR"

FD_RES = f"{FD_PRJ}/results"
FD_EXE = f"{FD_PRJ}/scripts"
FD_DAT = f"{FD_PRJ}/data"
FD_NBK = f"{FD_PRJ}/notebooks"
FD_DOC = f"{FD_PRJ}/docs"
FD_LOG = f"{FD_PRJ}/log"
FD_REF = f"{FD_PRJ}/references"

# =====================================
# helper function
# =====================================
def show_env():
    """ Quick print function to check variables"""
    print("BASE DIRECTORY (FD_BASE):", FD_BASE) 
    print("REPO DIRECTORY (FD_REPO):", FD_REPO)
    print("WORK DIRECTORY (FD_WORK):", FD_WORK)
    print("DATA DIRECTORY (FD_DATA):", FD_DATA)
    print("\n")
    
    print("You are working with     ", PROJECT) 
    print("PATH OF PROJECT (FD_PRJ):", FD_PRJ) 
    print("PROJECT RESULTS (FD_RES):", FD_RES) 
    print("PROJECT SCRIPTS (FD_EXE):", FD_EXE) 
    print("PROJECT DATA    (FD_DAT):", FD_DAT) 
    print("PROJECT NOTE    (FD_NBK):", FD_NBK) 
    print("PROJECT DOCS    (FD_DOC):", FD_DOC) 
    print("PROJECT LOG     (FD_LOG):", FD_LOG)  
    print("PROJECT REF     (FD_REF):", FD_REF) 
    print("")