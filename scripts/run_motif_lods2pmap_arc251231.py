"""
Batch precompute per-motif DP-based p-value mappings (PMAPs)
from .lods.pkl motif file.
"""

import numpy as np
import pickle
import os
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

from motifdelta import precompute_pmaps



# ====================================================================
# Function: loop through motifs and apply precompute_pmaps function
# --------------------------------------------------------------------
def batch_precompute_pmaps(
    dct_arr_motif_lod_WxB,
    arr_bg_B,
    num_alpha=1e-3,
    num_precision=1e-3,
    lst_score_range=(-50, 50),
    num_verbose_every=100
):
    """
    Batch precompute p-value mappings and T_bind thresholds for multiple motifs.

    Parameters
    ----------
    dct_arr_motif_lod_WxB : dict
        Dictionary of {motif_name: np.ndarray (W,B)} log-odds matrices.
        B = 4 if alphabet = A/C/G/T
    arr_bg_B : np.ndarray
        Background base probabilities (B,).
    num_alpha : float, optional
        Tail probability threshold for T_bind. Default = 1e-3.
    num_precision : float, optional
        Discretization bin width for DP grid. Default = 1e-3.
    lst_score_range : tuple, optional
        Range of scores to consider in bits. Default = (-50, 50).
    num_verbose_every : int, optional
        Print progress every N motifs.

    Returns
    -------
    dict
        {motif_name: precompute_pmaps(...) result dict}
    """
    ### init: total number of motifs and final results
    num_total = len(dct_arr_motif_lod_WxB)
    dct_results = {}

    for idx, (motif_name, arr_motif_lod_WxB) in enumerate(dct_arr_motif_lod_WxB.items()):
        try:
            dct_results[motif_name] = precompute_pmaps(
                arr_motif_lod_WxB,
                arr_bg_B,
                num_alpha=num_alpha,
                num_precision=num_precision,
                lst_score_range=lst_score_range
            )

        except Exception as e:
            print(f"[Warning] Skipped {motif_name} due to error: {e}")

        if (idx+1) % num_verbose_every == 0 or idx == num_total - 1:
            print(f"[{idx+1}/{num_total}] {motif_name}")

    print(f"Finished precomputing p-maps for {len(dct_results)} motifs.")
    return dct_results


def batch_precompute_pmaps_parallel(
    dct_arr_motif_lod_WxB,
    arr_bg_B,
    num_alpha=1e-3,
    num_precision=1e-3,
    lst_score_range=(-50, 50),
    num_verbose_every=100,
    num_workers=None
):
    """
    Parallel version of motif pmap precomputation.

    Parameters
    ----------
    
    Returns
    -------
    """
    ### init: set number of core for parallelization
    ### default: #{available CPU} - 1
    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 1)

    ### init: total number of motifs
    lst_motif_name = list(dct_arr_motif_lod_WxB.keys())
    num_total = len(lst_motif_name)

    ### init: final results with preserve order
    dct_results = {name: None for name in lst_motif_name}

    ### loop through each sequence and scan for all motifs in dictionary
    with ProcessPoolExecutor(max_workers=num_workers) as executor:

        ### assign each worker
        futures = {
            executor.submit(
                precompute_pmaps,
                arr_motif_lod_WxB,
                arr_bg_B,
                num_alpha,
                num_precision,
                lst_score_range,
            ): motif_name
            for motif_name, arr_motif_lod_WxB in dct_arr_motif_lod_WxB.items()
        }

        ### collect results
        for idx, future in enumerate(as_completed(futures)):
            ### get motif name
            motif_name = futures[future]

            ### get motif propbability mapping
            try:
                dct_results[motif_name] = future.result()
            ### error handling
            except Exception as e:
                print(f"[Warning] Skipped {motif_name} due to error: {e}")

            ### verbose
            if (idx+1) % num_verbose_every == 0 or idx == num_total - 1:
                print(f"[{idx+1}/{num_total}] {motif_name}")
                
    print(f"Finished precomputing p-maps for {len(dct_results)} motifs.")
    return dct_results


# ====================================================================
# Main function
# --------------------------------------------------------------------

def main(args):
    """Main function"""

    ### pass argument
    txt_fpath_inp = args.txt_fpath_inp
    txt_fpath_out = args.txt_fpath_out
    num_workers   = args.num_workers
    
    ### load motif
    with open(txt_fpath_inp, "rb") as f:
        obj = pickle.load(f)
    print(f"Loaded {txt_fpath_inp}")
    
    ### get motif log-odds and empirical background
    dct_arr_motif_lod_WxB = obj["lods"]
    arr_bg_B = obj["bg"]

    print(f"Loaded {len(dct_arr_motif_lod_WxB)} motifs")
    print(f"Loaded background: {arr_bg_B}")
    
    ### Run batch precomputation
    dct_results = batch_precompute_pmaps_parallel(
        dct_arr_motif_lod_WxB,
        arr_bg_B,
        num_alpha=1e-3,
        num_precision=1e-3,
        lst_score_range=(-40, 40),
        num_verbose_every=100,
        num_workers=num_workers
    )
    
    ### Save the dictionary
    with open(txt_fpath_out, "wb") as f:
        pickle.dump(dct_results, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Saved DP-maps for {len(dct_results)} motifs -> {txt_fpath_out}")

# ====================================================================
# CLI
# --------------------------------------------------------------------

if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Precompute DP-based p-value mappings (PMAPs) for motifs")
    
    parser.add_argument("--txt_fpath_inp",  type=str, required=True, help="Path to input motif pwd/lods file")
    parser.add_argument("--txt_fpath_out",  type=str, required=True, help="Path to output motif pmap file")
    parser.add_argument("--num_workers",    type=int, default =10,   help="Number of process for parallelization")
    
    args = parser.parse_args()

    ### run main function
    main(args)
