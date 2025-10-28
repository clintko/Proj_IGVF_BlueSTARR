"""
Batch precompute per-motif DP-based p-value mappings (PMAPs)
from JASPAR2024 .lods.npz motif file.
"""

import numpy as np
import os

from motifdelta import precompute_pmaps

# ====================================================================
# Function: loop through motifs and apply precompute_pmaps function
# --------------------------------------------------------------------
def batch_precompute_pmaps(
    dct_arr_motif_lod_Wx4,
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
    dct_arr_lod_WxB : dict
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
    dct_results = dict()
    num_total = len(dct_arr_lod_WxB)

    for idx, (motif_name, arr_motif_lod_WxB) in enumerate(dct_arr_motif_lod_WxB.items(), start=1):
        try:
            ### run recompute_pmaps
            dct_result = precompute_pmaps(
                arr_motif_lod_WxB,
                arr_bg_B,
                num_alpha=num_alpha,
                num_precision=num_precision,
                lst_score_range=lst_score_range
            )

            ### collect results
            dct_results[motif_name] = dct_result

            if num_verbose_every and (idx % num_verbose_every == 0 or idx == num_total):
                print(f"[{i}/{num_total}] Processed motif: {motif_name}")

        except Exception as e:
            print(f"[Warning] Skipped {motif_name} due to error: {e}")

    print(f"Finished precomputing p-maps for {len(dct_results)} motifs.")
    return dct_results


def batch_precompute_pmaps_parallel(
    dct_arr_motif_lod_Wx4,
    arr_bg_B,
    num_alpha=1e-3,
    num_precision=1e-3,
    lst_score_range=(-50, 50),
    num_verbose_every=100
):
    ### init: set number of core for parallelization
    ### default: #{available CPU} - 1
    num_total = len(dct_arr_motif_lod_WxB)
    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 1)


    ### init: final results
    dct_results = {}

    ### loop through each sequence and scan for all motifs in dictionary
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        ### assign workers
        futures = [
            executor.submit(scan_one_sequence, seq_record, dct_arr_motif_Wx4)
            for seq_record in lst_seq_record
        ]
        
        ### execute workers
        for f in as_completed(futures):
            idx, res = f.result()
            dct_results[idx] = res

    return dct_results

# ====================================================================
# Main function
# --------------------------------------------------------------------

def main(txt_fpath_inp, txt_fpath_out):
    """Main function"""

    ### load motif
    obj = np.load(txt_fpath_inp, allow_pickle=True)
    print(f"Loaded {txt_fpath_inp}")
    
    ### get motif log-odds and empirical background
    dct_arr_motif_lod_Wx4 = obj["lods"].item()
    arr_bg_B = obj["bg"]
    print(f"Loaded {len(dct_arr_lod_WxB)} motifs")
    print(f"Loaded background: {arr_bg_B}")
    
    ### Run batch precomputation
    dct_results = batch_precompute_pmaps(
        dct_arr_motif_lod_Wx4,
        arr_bg_B,
        num_alpha=1e-3,
        num_precision=1e-3,
        lst_score_range=(-40, 40),
        verbose_every=100
    )
    
    ### Save the dictionary
    np.savez_compressed(txt_fpath_out, dct_results=dct_results)
    print(f"Saved DP-maps for {len(dct_results)} motifs -> {txt_fpath_out}")
    
if __name__ == "__main__":

    ### Define input/output file path
    txt_fdiry_inp = "/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/results/analysis_variant_motif_richard"
    txt_fname_inp = "JASPAR2024_CORE_vertebrates_non-redundant.lods.npz"
    txt_fpath_inp = os.path.join(txt_fdiry_inp, txt_fname_inp)

    txt_fdiry_out = "/hpc/group/igvf/kk319/repo/Proj_IGVF_BlueSTARR/results/analysis_variant_motif_richard"
    txt_fname_out = "JASPAR2024_CORE_vertebrates_non-redundant.pmap.npz"
    txt_fpath_out = os.path.join(txt_fdiry_out, txt_fname_out)
    
    ### Run main function
    main(txt_fpath_inp, txt_fpath_out)