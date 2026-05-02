"""
Import motif PWM pickle and convert PWMs into log-odds scores
"""

import numpy as np
import pickle
import argparse

from motifdelta.score import pwm_to_logodds


# ====================================================================
# Main function
# --------------------------------------------------------------------

def main(args):
    """Main function"""

    ### pass argument
    txt_fpath_inp = args.txt_fpath_inp
    txt_fpath_bkg = args.txt_fpath_bkg
    txt_fpath_out = args.txt_fpath_out
        
    ### Load empirical background
    arr_bkg_B = np.load(txt_fpath_bkg)
    arr_bkg_B = np.asarray(arr_bkg_B, dtype=float).reshape(-1)
    arr_bkg_B = arr_bkg_B / arr_bkg_B.sum()
    print(f"Loaded background: {arr_bkg_B}")

    ### load input motif pwm
    with open(txt_fpath_inp, "rb") as file:
        dct_arr_motif_pwm_WxB = pickle.load(file)
    print(f"Loaded {len(dct_arr_motif_pwm_WxB)} motifs")
    print("Example keys:", list(dct_arr_motif_pwm_WxB)[:3])
    
    first_key = next(iter(dct_arr_motif_pwm_WxB))
    first_val = dct_arr_motif_pwm_WxB[first_key]
    print("First PWM shape:", np.asarray(first_val).shape)
    
    ### init: log-odds
    dct_arr_motif_lod_WxB = dict()

    for txt_motif_name, arr_motif_pwm_WxB in dct_arr_motif_pwm_WxB.items():
        ### Convert pwm (W,4) -> log-odds (W,4)
        arr_motif_lod_WxB = pwm_to_logodds(arr_motif_pwm_WxB, arr_bkg_B)

        ### Collect results
        dct_arr_motif_lod_WxB[txt_motif_name] = arr_motif_lod_WxB
    
    ### save results
    dct_out = {
        "pwms": dct_arr_motif_pwm_WxB,
        "lods": dct_arr_motif_lod_WxB,
        "bg":   arr_bkg_B,
        "names": list(dct_arr_motif_pwm_WxB.keys()),
        "alphabet": "ACGT"
    }

    print(f"Saved {len(dct_arr_motif_pwm_WxB)} motifs -> {txt_fpath_out}")    
    with open(txt_fpath_out, "wb") as file:
        pickle.dump(dct_out, file, protocol=pickle.HIGHEST_PROTOCOL)

# ====================================================================
# CLI
# --------------------------------------------------------------------

if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files")
    
    parser.add_argument("--txt_fpath_inp",  type=str, required=True, help="Path to input  motif pwm pickle file")
    parser.add_argument("--txt_fpath_bkg",  type=str, required=True, help="Path to input  motif background file")
    parser.add_argument("--txt_fpath_out",  type=str, required=True, help="Path to output motif lod pickle file")
    
    args = parser.parse_args()

    ### run main function
    main(args)

    