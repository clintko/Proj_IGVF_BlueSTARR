"""
Module for DP-based p-value mapping and T_bind calculation
Note: T_bind is the critical value for given a significance level alpha
"""

import numpy as np


def build_score_distribution(arr_lod_WxB, arr_bg_B, num_precision=1e-2, lst_score_range=(-50, 50)):
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
    lst_score_range : tuple, optional
        Range of scores to consider (min, max) in bits.

    Returns
    -------
    arr_num_score_grid : np.ndarray
        Discretized score grid (bins) in bits.
    arr_num_pmf : np.ndarray
        Probability mass function P(S = s) evaluated on 'arr_num_score_grid'.
    """
    ### sanity check
    if arr_bg_B.shape[0] != arr_lod_WxB.shape[1]:
        raise ValueError("Background length B does not match motif columns.")
        
    if not np.isclose(arr_bg_B.sum(), 1.0):
        raise ValueError("Background probabilities must sum to 1.")
    
    ### initialize a score grid
    num_step = num_precision
    arr_num_score_grid = np.arange(
        lst_score_range[0], lst_score_range[1] + num_step, num_step
    )
    num_bins = len(arr_num_score_grid)

    ### start with a degenerate (Dirac delta) distribution at score = 0
    arr_num_score_pmf = np.zeros(num_bins)
    idx_zero          = np.searchsorted(arr_num_score_grid, 0)
    
    if idx_zero < num_bins:
        arr_num_score_pmf[idx_zero] = 1.0
    else:
        raise ValueError("Score range too small; zero not included.")

    ### dynamic programming convolution across motif positions
    for idx_pos in range(arr_lod_WxB.shape[0]):
        ### init a new pmf
        arr_num_new_pmf = np.zeros_like(arr_num_score_pmf)
        
        for idx_base in range(arr_lod_WxB.shape[1]):  # usually 4 bases
            
            ### convolution to get the update pmf
            num_shift = int(round(arr_lod_WxB[idx_pos, idx_base] / num_step))
            #arr_num_new_pmf += arr_bg_B[idx_base] * np.roll(arr_num_score_pmf, num_shift)
            
            w = arr_bg_B[idx_base]
            
            ### skip if shift moves everything out of range
            if abs(num_shift) >= num_bins:    
                continue
            if num_shift > 0:
                arr_num_new_pmf[num_shift:] += w * arr_num_score_pmf[:-num_shift]
            elif num_shift < 0:
                s = -num_shift
                arr_num_new_pmf[:-s] += w * arr_num_score_pmf[s:]
            else:
                arr_num_new_pmf += w * arr_num_score_pmf
        
        ### normalize into a valid pmf
        arr_num_score_pmf = arr_num_new_pmf
        num_score_sum = arr_num_score_pmf.sum()
        if num_score_sum <= 0:
            raise RuntimeError("PMF became degenerate (sum of pmf == 0); consider widen score range or reduce precision.")
        arr_num_score_pmf /= num_score_sum

    return arr_num_score_grid, arr_num_score_pmf


def build_score_to_pvalue(arr_num_score_grid, arr_num_score_pmf):
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
    arr_num_score_ccdf = np.cumsum(arr_num_score_pmf[::-1])[::-1]

    ### normalized the max into as probability one
    arr_num_score_ccdf /= arr_num_score_ccdf[0]  
    return arr_num_score_grid, arr_num_score_ccdf


def find_Tbind(arr_num_score_grid, arr_num_score_ccdf, num_alpha=1e-3):
    """
    Find the score threshold T_bind such that P(S >= T_bind) = alpha.

    Parameters
    ----------
    arr_num_score_grid : np.ndarray
        Discretized score grid.
    arr_num_score_ccdf : np.ndarray
        complementary cumulative distribution function from build_score_to_pvalue().
    num_alpha : float
        Desired tail probability threshold.

    Returns
    -------
    float
        Score threshold corresponding to P(S >= s) ~ alpha.
    """
    ### find the smallest score index where ccdf = 1 - CDF <= alpha
    ### i.e. critical values of corresponding PMF
    
    ### reverse since CCDF is descending
    arr_num_score_ccdf_rev = arr_num_score_ccdf[::-1]

    ### search through the acending array
    idx_rev = np.searchsorted(arr_num_score_ccdf_rev, num_alpha, side="left")

    ### get the exact index of the critical value/threshold
    idx = len(arr_num_score_ccdf) - idx_rev - 1

    ### clamp index and return threshold
    idx = np.clip(idx, 0, len(arr_num_score_grid) - 1)
    return float(arr_num_score_grid[idx])


def precompute_pmaps(arr_lod_WxB, arr_bg_B, num_alpha=1e-3, num_precision=1e-2, lst_score_range=(-60, 60)):
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
            "arr_num_score_pmf":  np.ndarray,  # PMF
            "arr_num_score_ccdf": np.ndarray,  # right-tail probability (P(S >= s))
            "arr_lod_WxB": np.ndarray,         # the motif lod used to generate teh probability map
            "num_Tbind": float,                # score threshold for alpha (critical value)
            "num_alpha": float,                # alpha used for thresholding
            "num_precision": float             # bin width
        }
    """
    ### build PMF
    arr_num_score_grid, arr_num_score_pmf = build_score_distribution(
        arr_lod_WxB,
        arr_bg_B,
        num_precision   = num_precision,
        lst_score_range = lst_score_range
    )

    ### convert to right-tail CCDF for p-value mapping
    arr_num_score_grid, arr_num_score_ccdf = build_score_to_pvalue(
        arr_num_score_grid,
        arr_num_score_pmf
    )

    ### find threshold T_bind for desired alpha using Complement CDF (CCDF)
    num_Tbind = find_Tbind(
        arr_num_score_grid,
        arr_num_score_ccdf,
        num_alpha=num_alpha
    )

    ### collect results 
    dct_result = {
        "arr_num_score_grid": arr_num_score_grid,
        "arr_num_score_pmf":  arr_num_score_pmf,
        "arr_num_score_ccdf": arr_num_score_ccdf,
        #"arr_lod_WxB":   arr_lod_WxB,
        "num_Tbind":     num_Tbind,
        "num_alpha":     num_alpha,
        "num_precision": num_precision
    }

    return dct_result

def map_score_to_pvalue(
    arr_num_score,
    arr_num_score_grid,
    arr_num_score_ccdf,
    do_interpolate = False,
    do_clip = True,
):
    """
    Map motif scores to right-tail p-values using a precomputed score grid and CCDF.

    Parameters
    ----------
    arr_num_score : np.ndarray
        Vector/array of motif scores to map (e.g., from scan_motif / scan_motif_both_strands).
    arr_num_score_grid : np.ndarray
        Discretized score grid (strictly increasing) used to build the CCDF.
    arr_num_score_ccdf : np.ndarray
        Right-tail probabilities P(S >= s) aligned to arr_num_score_grid.
    do_interpolate : bool, optional
        If True, linearly interpolate p-values between grid points.
        If False, use left-binning via searchsorted (conservative).
    do_clip : bool, optional
        If True, clip scores to [grid_min, grid_max] before mapping.
        If False, scores outside the grid will be snapped to nearest edge.

    Returns
    -------
    np.ndarray
        Right-tail p-values corresponding to arr_num_score.
    """
    ### Optional clip (avoid index warnings)
    if do_clip:
        smin, smax = arr_num_score_grid[0], arr_num_score_grid[-1]
        arr_score = np.clip(arr_num_score, smin, smax)
    else:
        arr_score = arr_num_score

    if not do_interpolate:
        ### Step mapping via searchsorted (left): P(S >= s) ~ CCDF[bin(s)]
        idx = np.searchsorted(arr_num_score_grid, arr_score, side="left")
        idx = np.clip(idx, 0, len(arr_num_score_grid) - 1)
        arr_num_pval = arr_num_score_ccdf[idx]
        return arr_num_pval

    ### Linear interpolation: smooths p-values between neighboring bins
    ### np.interp assumes x is increasing (grid is increasing)
    arr_num_pval = np.interp(arr_score, arr_num_score_grid, arr_num_score_ccdf)
    return arr_num_pval