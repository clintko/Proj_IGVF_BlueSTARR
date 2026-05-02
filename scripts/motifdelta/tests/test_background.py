"""
Unit tests for background model calculation
"""

import numpy as np
from motifdelta.stats import background_zero_order

def test_background_case01_uniform():
    """Uniform background across A, C, G, T"""
    seqs = ["ACGTACGT", "AAAA", "CCCC", "GGGG", "TTTT"]
    bg = background_zero_order(seqs)
    assert np.isclose(bg.sum(), 1.0, atol=1e-8)
    assert np.allclose(bg, 0.25, atol=0.01)


def test_background_case02_at_rich():
    """AT-rich background example"""
    seqs = ["AAAAATTTTT", "ATATATATAT", "AAAATTTT"]
    bg = background_zero_order(seqs)
    gc = bg[1] + bg[2]  # C + G; default order: ACGT
    assert np.isclose(bg.sum(), 1.0, atol=1e-8)
    assert gc < 0.5  # should be AT-rich