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
#import sys
#from pathlib import Path

#FD_EXE = Path("/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/scripts")
#if str(FD_EXE) not in sys.path:
#    sys.path.insert(0, str(FD_EXE))

from motifdelta import pwm_to_logodds


# ====================================================================
# Main function
# --------------------------------------------------------------------

def main(args):
    """Main function"""

    ### pass argument
    txt_fpath_inp = args.txt_fpath_inp
    txt_fpath_bgk = args.txt_fpath_bgk
    txt_fpath_out = args.txt_fpath_out
        
    ### Load empirical background
    arr_bg_B = np.load(txt_fpath_bgk)
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
    
    ### save results
    dct_out = {
        "pwms": dct_arr_motif_pwm_Wx4,
        "lods": dct_arr_motif_lod_Wx4,
        "bg": arr_bg_B,
        "names": list(dct_arr_motif_pwm_Wx4.keys()),
        "alphabet": "ACGT"
    }

    print(f"Saved {len(dct_arr_motif_pwm_Wx4)} motifs -> {txt_fpath_out}")    
    with open(txt_fpath_out, "wb") as f:
        pickle.dump(dct_out, f, protocol=pickle.HIGHEST_PROTOCOL)


# ====================================================================
# CLI
# --------------------------------------------------------------------

if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files")
    
    parser.add_argument("--txt_fpath_inp",  type=str, required=True, help="Path to input motif meme file")
    parser.add_argument("--txt_fpath_bgk",  type=str, required=True, help="Path to input motif background file")
    parser.add_argument("--txt_fpath_out",  type=str, required=True, help="Path to output motif pickle file")
    
    args = parser.parse_args()

    ### run main function
    main(args)

    