"""
Unit tests for deltamotif.model (DP-based score distribution & T_bind)
"""

import numpy as np
from motifdelta.model import (
    build_score_distribution,
    build_score_to_pvalue,
    find_Tbind,
    precompute_pmaps
)


def test_build_score_distribution_sum_to_one():
    """PMF should always sum to 1 (normalized probability)."""
    arr_lod_WxB = np.array([[1, -1, 0, 0]])  # simple 1-position motif
    arr_bg_B = np.array([0.25, 0.25, 0.25, 0.25])
    arr_score, arr_pmf = build_score_distribution(arr_lod_WxB, arr_bg_B, num_precision=0.1)
    assert np.isclose(arr_pmf.sum(), 1.0, atol=1e-6)


def test_build_score_to_pvalue_monotonic():
    """CCDF (right-tail) should be monotonically non-increasing."""
    arr_lod_WxB = np.array([[1, 0, 0, -1]])
    arr_bg_B = np.array([0.25, 0.25, 0.25, 0.25])
    arr_score, arr_pmf = build_score_distribution(arr_lod_WxB, arr_bg_B, num_precision=0.1)
    arr_score, arr_ccdf = build_score_to_pvalue(arr_score, arr_pmf)
    assert np.all(np.diff(arr_ccdf) <= 1e-9)  # monotonic decreasing


def test_find_Tbind_threshold_moves_with_alpha():
    """T_bind should increase (stricter) when alpha decreases."""
    arr_lod_WxB = np.array([[1.0, -0.5, 0.2, -1.0], [-0.5, 0.5, 0.8, -0.2]])
    arr_bg_B = np.array([0.25, 0.25, 0.25, 0.25])
    arr_score, arr_pmf = build_score_distribution(arr_lod_WxB, arr_bg_B, num_precision=0.05)
    arr_score, arr_ccdf = build_score_to_pvalue(arr_score, arr_pmf)

    T1 = find_Tbind(arr_score, arr_ccdf, num_alpha=1e-2)
    T2 = find_Tbind(arr_score, arr_ccdf, num_alpha=1e-4)
    assert T2 >= T1  # smaller alpha -> same or higher threshold


def test_precompute_pmaps_output_keys():
    """Precompute dict should include all expected keys."""
    arr_lod_WxB = np.array([[1.0, -1.0, 0.5, -0.5], [-0.5, 0.5, 1.0, -1.0]])
    arr_bg_B = np.array([0.25, 0.25, 0.25, 0.25])
    dct = precompute_pmaps(arr_lod_WxB, arr_bg_B, num_alpha=1e-3, num_precision=0.1)

    expected_keys = {
        "arr_num_score_grid",
        "arr_num_score_pmf",
        "arr_num_score_ccdf",
        "arr_lod_WxB",
        "num_Tbind",
        "num_alpha",
        "num_precision"
    }
    assert set(dct.keys()) == expected_keys
    assert np.isclose(dct["arr_num_score_pmf"].sum(), 1.0, atol=1e-6)
    assert np.all(np.diff(dct["arr_num_score_ccdf"]) <= 1e-9)
    assert isinstance(dct["num_Tbind"], float)


def test_precompute_pmaps_Tbind_reasonable_range():
    """T_bind should lie within the score grid range."""
    arr_lod_WxB = np.array([[1, -1, 1, -1]])
    arr_bg_B = np.array([0.3, 0.2, 0.2, 0.3])
    dct = precompute_pmaps(arr_lod_WxB, arr_bg_B)
    min_score, max_score = dct["arr_num_score_grid"][0], dct["arr_num_score_grid"][-1]
    Tbind = dct["num_Tbind"]
    assert min_score <= Tbind <= max_score
