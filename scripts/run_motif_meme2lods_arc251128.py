"""
Import meme file and convert pwm into log-odds scores
"""

import numpy as np
import pickle
import time
import os
import argparse

from memelite.io import read_meme

# ====================================================================
# Import pwm_to_logodds from motifdelta
# --------------------------------------------------------------------
import sys
from pathlib import Path

FD_EXE = Path("/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/scripts")
if str(FD_EXE) not in sys.path:
    sys.path.insert(0, str(FD_EXE))

from motifdelta import pwm_to_logodds

# ====================================================================
# Main function
# --------------------------------------------------------------------

def main(txt_fpath_inp, txt_fpath_out, txt_fpath_bg):
    """Main function"""

    ### Load empirical background
    arr_bg_B = np.load(txt_fpath_bg)
    print(f"Loaded background: {arr_bg_B}")

    ### Read MEME motifs
    dct_arr_motif_pwm_4xW = read_meme(txt_fpath_inp)
    print(f"Loaded {len(dct_arr_motif_pwm_4xW)} motifs from MEME file")

    ### init: PWM and log-odds
    dct_arr_motif_pwm_Wx4 = dict()
    dct_arr_motif_lod_Wx4 = dict()

    for motif_name, arr_motif_pwm_4xW in dct_arr_motif_pwm_4xW.items():
        ### Convert pwm (4,W) -> pwm (W,4) -> log-odds (W,4)
        arr_pwm_Wx4 = arr_motif_pwm_4xW.T.astype(float)
        arr_lod_Wx4 = pwm_to_logodds(arr_pwm_Wx4, arr_bg_B)

        ### Collect results
        dct_arr_motif_pwm_Wx4[motif_name] = arr_pwm_Wx4
        dct_arr_motif_lod_Wx4[motif_name] = arr_lod_Wx4
    
    ### Export motif PWM and log-odds
    #np.savez_compressed(
    #    txt_fpath_out,
    #    pwms  = np.array(dct_arr_motif_pwm_Wx4, dtype=object),
    #    lods  = np.array(dct_arr_motif_lod_Wx4, dtype=object),
    #    bg    = np.array(arr_bg_B, dtype=np.float32),
    #    names = np.array(list(dct_arr_motif_pwm_Wx4.keys()), dtype=object),
    #    alphabet="ACGT",
    #)

    dct_save = {
        "pwms": dct_arr_motif_pwm_Wx4,
        "lods": dct_arr_motif_lod_Wx4,
        "bg": arr_bg_B,
        "names": list(dct_arr_motif_pwm_Wx4.keys()),
        "alphabet": "ACGT"
    }

    print(f"Saved {len(dct_arr_motif_pwm_Wx4)} motifs -> {txt_fpath_out}")    
    with open(txt_fpath_out, "wb") as f:
        pickle.dump(dct_save, f, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":

    ### Define input/output file path
    txt_fdiry_inp = "/hpc/group/igvf/kk319/data/jaspar2024"
    txt_fname_inp = "JASPAR2024_CORE_vertebrates_non-redundant.meme"
    txt_fpath_inp = os.path.join(txt_fdiry_inp, txt_fname_inp)

    txt_fdiry_out = "/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/results/analysis_variant_motif_richard"
    txt_fname_out = "JASPAR2024_CORE_vertebrates_non-redundant.lods.pkl"
    txt_fpath_out = os.path.join(txt_fdiry_out, txt_fname_out)

    txt_fname_bg  = "background_zero_order.npy"
    txt_fpath_bg  = os.path.join(txt_fdiry_out, txt_fname_bg)
    
    ### Run main function
    main(txt_fpath_inp, txt_fpath_out, txt_fpath_bg)
    