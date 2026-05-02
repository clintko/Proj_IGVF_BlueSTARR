"""
Module for scanning DNA sequences with motif matrices
"""

import numpy as np
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

from motifdelta.encode import one_hot_encode
from motifdelta.rc     import reverse_complement_matrix
from motifdelta.model  import map_score_to_pvalue


def scan_motif(txt_seq, arr_motif_WxB, txt_alphabet="ACGT"):
    """
    Scan a motif across a DNA sequence using sliding windows.

    Parameters
    ----------
    txt_seq : str
        DNA sequence (string of A/C/G/T) with length L
    arr_motif_WxB : np.ndarray
        Motif score matrix of shape (W,B) (e.g. motif pwm or motif log-odds)
    txt_alphabet : str
        Alphabet order for both sequence and motif (default: ACGT)

    Returns
    -------
    np.ndarray
        Vector of motif scores with shape (L - W + 1,)
    """
    ### encode sequence
    X = one_hot_encode(txt_seq, txt_alphabet=txt_alphabet)
    L = X.shape[0]
    W = arr_motif_WxB.shape[0]

    if L < W:
        raise ValueError(f"Sequence length {L} shorter than motif length {W}")

    ### sliding window view: shape (L-W+1, W, B)
    arr_windows = np.lib.stride_tricks.sliding_window_view(X, window_shape=(W, X.shape[1]), axis=(0,1))
    arr_windows = arr_windows.reshape(-1, W, X.shape[1])
    
    ### element-wise multiply and sum over both dimensions
    arr_scores = np.tensordot(arr_windows, arr_motif_WxB, axes=([1, 2], [0, 1]))

    return arr_scores

def scan_motif_both_strands(txt_seq, arr_motif_WxB, txt_alphabet="ACGT"):
    """
    Scan motif on both forward and reverse-complement strands.

    Returns
    -------
    tuple (arr_score_fwd, arr_score_rev)
        arr_score_fwd : np.ndarray
            Motif scores on forward strand
        arr_score_rev : np.ndarray
            Motif scores on reverse-complement strand
    """
    ### forward scanning
    arr_score_fwd = scan_motif(txt_seq, arr_motif_WxB, txt_alphabet=txt_alphabet)

    ### reverse scanning
    arr_motif_rc  = reverse_complement_matrix(arr_motif_WxB, txt_alphabet=txt_alphabet)
    arr_score_rev = scan_motif(txt_seq, arr_motif_rc, txt_alphabet=txt_alphabet)
    
    return arr_score_fwd, arr_score_rev

def scan_one_sequence(
    txt_seq_idx: str,
    txt_seq_str: str,
    dct_motif_model,
    txt_alphabet="ACGT"
):
    """
    Scan all motifs on a single sequence (both strands), returning scores, p-values
    This function is for a single sequence only (no delta here).

    Parameters
    ----------
    txt_seq_idx: str,
    txt_seq_str: str,
    dct_motif_model: dict
            "arr_num_score_grid": np.ndarray,  # score bins
            "arr_num_score_pmf":  np.ndarray,  # score PMF
            "arr_num_score_ccdf": np.ndarray,  # score right-tail probability (P(S >= s))
            "num_Tbind": float,                # score threshold for alph (critical value)
            "num_alpha": float,                # alpha used for thresholding
            "num_precision": float             # bin width
    txt_alphabet: str
        Default = "ACGT"
        
    Returns
    -------
    
    """
    ### init: collect results
    dct_results = {}

    ### loop through each motif
    for txt_motif_name, dct_motif_pmap in dct_motif_model.items():
        ### 
        arr_lod_WxB        = dct_motif_pmap["arr_lod_WxB"] 
        arr_num_score_grid = dct_motif_pmap["arr_num_score_grid"]
        arr_num_score_ccdf = dct_motif_pmap["arr_num_score_ccdf"]
        num_Tbind          = dct_motif_pmap["num_Tbind"]

        ###
        arr_num_score_fwd, arr_num_score_rev = scan_motif_both_strands(
            txt_seq_str, arr_lod_WxB, txt_alphabet=txt_alphabet
        )

        ### vectorized mapping score to p-values using searchsorted
        #idx_fwd = np.searchsorted(arr_num_score_grid, arr_num_score_fwd, side="left")
        #idx_rev = np.searchsorted(arr_num_score_grid, arr_num_score_rev, side="left")
        #idx_fwd = np.clip(idx_fwd, 0, len(arr_num_score_grid)-1)
        #idx_rev = np.clip(idx_rev, 0, len(arr_num_score_grid)-1)
        #arr_num_pval_fwd = arr_num_score_ccdf[idx_fwd]
        #arr_num_pval_rev = arr_num_score_ccdf[idx_rev]

        ### forward strand
        arr_num_pval_fwd = map_score_to_pvalue(
            arr_num_score_fwd,
            arr_num_score_grid,
            arr_num_score_ccdf,
            do_interpolate=True,
            do_clip=True
        )
        
        ### reverse strand
        arr_num_pval_rev = map_score_to_pvalue(
            arr_num_score_rev,
            arr_num_score_grid,
            arr_num_score_ccdf,
            do_interpolate=True,
            do_clip=True
        )
        
        ### loop through each motif and scan
        dct_results[txt_motif_name] = {
            "arr_num_score_forward": arr_num_score_fwd,
            "arr_num_score_reverse": arr_num_score_rev,
            "arr_num_pval_forward":  arr_num_pval_fwd,
            "arr_num_pval_reverse":  arr_num_pval_rev,
            "num_Tbind":             num_Tbind
        }

    return txt_seq_idx, dct_results

def scan_sequence_batch_parallel(
    lst_txt_seq_idx,
    lst_txt_seq_str,
    dct_motif_model,
    num_workers=None
):
    """
    Parallel scan using precomputed motif p-maps (scores + p-values).

    Parameters
    ----------
    lst_txt_seq_idx : list[str]
        List of sequence IDs.
    lst_txt_seq_str : list[str]
        List of DNA sequences (ACGT), same length/order as lst_txt_seq_idx.
    dct_motif_model : dict
        motif_name -> {
            "arr_lod_WxB", "arr_num_score_grid", "arr_num_score_ccdf", "num_Tbind", ...
        }
    num_workers : int, optional
        Number of worker processes. Defaults to CPU count - 1.
    """
    if len(lst_txt_seq_idx) != len(lst_txt_seq_str):
        raise ValueError("lst_txt_seq_idx and lst_txt_seq_str must have the same length.")

    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 1)

    dct_results = {}
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(
                scan_one_sequence,
                txt_seq_idx,
                txt_seq_str,
                dct_motif_model
            )
            for txt_seq_idx, txt_seq_str in zip(lst_txt_seq_idx, lst_txt_seq_str)
        ]

        for f in as_completed(futures):
            txt_seq_idx, res = f.result()
            dct_results[txt_seq_idx] = res

    return dct_results

