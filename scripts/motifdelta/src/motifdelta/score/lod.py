"""
Module for motif score calculation
"""

import numpy as np

def pwm_to_logodds(arr_pwm_WxB, arr_bkg_B=(0.25, 0.25, 0.25, 0.25), eps=1e-6):
    """
    Convert a motif PWM (W,B) to a log-likelihood ratio / log-odds scoring matrix (W,B)

    Parameters
    ----------
    arr_pwm_WxB : np.ndarray
        Motif PWM matrix of shape (W,B), where rows are positions
        and columns correspond to bases (A,C,G,T).
    arr_bkg_B : array-like, shape (B,)
        Background base probabilities (default: uniform [0.25,0.25,0.25,0.25]).
    eps : float
        Small constant added to avoid log(0).

    Returns
    -------
    np.ndarray
        Motif log-odds matrix of shape (W,B), in bits.
    """
    ### clip to avoid zeros
    arr_pwm_WxB = np.clip(np.asarray(arr_pwm_WxB, dtype=float), eps, 1.0)
    arr_bkg_B   = np.clip(np.asarray(arr_bkg_B,    dtype=float), eps, 1.0)

    ### normalize background to ensure it sums to 1
    arr_bkg_B /= arr_bkg_B.sum()

    ### calculate log-odds score
    arr_lod_WxB = np.log2(arr_pwm_WxB / arr_bkg_B[None, :])
    return arr_lod_WxB