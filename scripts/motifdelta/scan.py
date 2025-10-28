"""
Module for scanning DNA sequences with motif matrices
"""

import numpy as np
from motifdelta.encode import one_hot_encode
from motifdelta.rc     import reverse_complement_matrix


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
    arr_windows = np.lib.stride_tricks.sliding_window_view(X, (W, X.shape[1]))
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

