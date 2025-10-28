"""
Module for DP-based p-value mapping and T_bind calculation
"""

import numpy as np


def build_score_distribution(arr_lod_WxB, arr_bg_B, num_precision=1e-3, lst_score_range=(-50, 50)):
    """
    Compute discrete probability mass function (PMF) of motif scores under
    a 0th-order background model using dynamic programming.

    Parameters
    ----------
    arr_lod_WxB : np.ndarray
        Motif log-odds matrix of shape (W,4).
    arr_bg_B : np.ndarray
        Background base probabilities, shape (4,).
    num_precision : float, optional
        Discretization step size (bin width in bits). Default = 1e-3.
    tup_score_range : tuple, optional
        Range of scores to consider (min, max) in bits.

    Returns
    -------
    arr_num_score_grid : np.ndarray
        Discretized score grid (bins) in bits.
    arr_num_pmf : np.ndarray
        Probability mass function P(S = s) evaluated on 'arr_num_score_grid'.
    """
    ### initialize a score grid
    num_step = num_precision
    arr_num_score_grid = np.arange(
        lst_score_range[0], lst_score_range[1] + num_step, num_step
    )
    num_bins = len(arr_num_score_grid)

    ### start with a degenerate (Dirac delta) distribution at score = 0
    arr_num_pmf = np.zeros(num_bins)
    idx_zero    = np.searchsorted(arr_num_score_grid, 0)
    
    if idx_zero < num_bins:
        arr_num_pmf[idx_zero] = 1.0
    else:
        raise ValueError("Score range too small; zero not included.")

    ### dynamic programming convolution across motif positions
    for idx_pos in range(arr_lod_WxB.shape[0]):
        ### init a new pmf
        arr_num_new_pmf = np.zeros_like(arr_num_pmf)
        
        for idx_base in range(arr_lod_WxB.shape[1]):  # usually 4 bases
            ### convolution to get the update pmf
            num_shift = int(round(arr_lod_WxB[idx_pos, idx_base] / num_step))
            arr_num_new_pmf += arr_bg_B[idx_base] * np.roll(arr_num_pmf, num_shift)
        
        ### normalize into a valid pmf
        arr_num_pmf = arr_num_new_pmf
        arr_num_pmf /= arr_num_pmf.sum()

    return arr_num_score_grid, arr_num_pmf



def build_score_to_pvalue(arr_num_score_grid, arr_num_pmf):
    """
    Convert PMF to a right-tail cdf function (complementary cumulative distribution function; for score to p-value mapping).

    Parameters
    ----------
    scores : np.ndarray
        Score grid from build_score_distribution().
    pmf : np.ndarray
        Corresponding probability mass function.

    Returns
    -------
    scores : np.ndarray
        Same as input.
    ccdf : np.ndarray
        P(S >= s) for each score.
    """
    ### cumulative distribution from the left
    arr_num_ccdf = np.cumsum(arr_num_pmf[::-1])[::-1]

    ### normalized the max into as probability one
    arr_num_ccdf /= arr_num_ccdf[0]  
    return arr_num_score_grid, arr_num_ccdf


def find_Tbind(arr_num_score_grid, arr_num_ccdf, num_alpha=1e-3):
    """
    Find the score threshold T_bind such that P(S >= T_bind) = alpha.

    Parameters
    ----------
    scores : np.ndarray
        Discretized score grid.
    sf : np.ndarray
        complementary cumulative distribution function from build_score_to_pvalue().
    alpha : float
        Desired tail probability threshold.

    Returns
    -------
    float
        Score threshold corresponding to P(S >= s) ~ alpha.
    """
    ### find the smallest score index where ccdf = 1 - CDF <= alpha
    ### i.e. critical values of corresponding PMF
    
    ### previous code
    #idx = np.searchsorted(arr_num_ccdf, num_alpha, side="left")
    #idx = min(idx, len(arr_num_score_grid) - 1)
    #return arr_num_score_grid[idx]

    ### reverse since CCDF is descending
    arr_ccdf_rev = arr_num_ccdf[::-1]

    ### search through the acending array
    idx_rev = np.searchsorted(arr_ccdf_rev, num_alpha, side="left")

    ### get the exact index of the critical value/threshod
    idx = len(arr_num_ccdf) - idx_rev - 1

    ### clamp index and return threshold
    idx = np.clip(idx, 0, len(arr_num_score_grid) - 1)
    return arr_num_score_grid[idx]


def precompute_pmaps(arr_lod_WxB, arr_bg_B, num_alpha=1e-3, num_precision=1e-3, lst_score_range=(-50, 50)):
    """
    Precompute p-value mapping and T_bind threshold for a single motif
    based on its log-odds matrix and background model.

    Parameters
    ----------
    arr_lod_WxB : np.ndarray
        Motif log-odds matrix (W,B). B is 4 if alphabet = A/C/G/T
    arr_bg_B : np.ndarray
        Background base probabilities (B,).
    num_alpha : float, optional
        Desired right-tail probability threshold. Default = 1e-3.
    num_precision : float, optional
        Discretization bin width for score distribution. Default = 1e-3.
    lst_score_range : list or tuple, optional
        Range of scores to consider (min, max). Default = (-50, 50).

    Returns
    -------
    dct_result : dict
        {
            "arr_num_score_grid": np.ndarray,  # score bins
            "arr_num_pmf":  np.ndarray,        # PMF
            "arr_num_ccdf": np.ndarray,        # right-tail probability (P(S >= s))
            "num_Tbind": float,                # score threshold for alph (critical value)
            "num_alpha": float,                # alpha used for thresholding
            "num_precision": float             # bin width
        }
    """
    ### build PMF
    arr_num_score_grid, arr_num_pmf = build_score_distribution(
        arr_lod_WxB,
        arr_bg_B,
        num_precision   = num_precision,
        lst_score_range = lst_score_range
    )

    ### convert to right-tail CCDF for p-value mapping
    arr_num_score_grid, arr_num_ccdf = build_score_to_pvalue(
        arr_num_score_grid,
        arr_num_pmf
    )

    ### find threshold T_bind for desired alpha using Complement CDF (CCDF)
    num_Tbind = find_Tbind(
        arr_num_score_grid,
        arr_num_ccdf,
        num_alpha=num_alpha
    )

    ### collect results 
    dct_result = {
        "arr_num_score_grid": arr_num_score_grid,
        "arr_num_pmf":  arr_num_pmf,
        "arr_num_ccdf": arr_num_ccdf,
        "num_Tbind": num_Tbind,
        "num_alpha": num_alpha,
        "num_precision": num_precision
    }

    return dct_result
