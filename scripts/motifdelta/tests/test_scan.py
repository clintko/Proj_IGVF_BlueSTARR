"""
Unit tests for fast batched motif scanning functions (motifdelta.score.scan)
"""

import numpy as np
import pytest

from motifdelta.seq   import reverse_complement_sequence, reverse_complement_matrix
from motifdelta.score import (
    prepare_motif_kernels,
    one_hot_batch,
    scan_sequence_batch,
)


def _scan_single_sequence(txt_seq: str, arr_motif_WxB: np.ndarray, txt_alphabet: str = "ACGT"):
    """
    Helper: scan a single sequence with a single motif using the batched API.
    Returns forward/reverse score vectors of shape (P,).
    """
    d = {"m": np.asarray(arr_motif_WxB, dtype=np.float32)}
    arr_fwd, arr_rev, names, W_max = prepare_motif_kernels(d, txt_alphabet=txt_alphabet)

    X = one_hot_batch([txt_seq], txt_alphabet=txt_alphabet)  # (1, L, B)
    out = scan_sequence_batch(X, arr_fwd, arr_rev)           # (1, 1, P, 2)

    return out[0, 0, :, 0], out[0, 0, :, 1]


def test_scan_basic_shape_and_finite():
    txt_seq = "ACGTAC"
    arr_motif = np.array([
        [ 1, -1, -1, -1],  # A
        [-1,  1, -1, -1],  # C
        [-1, -1,  1, -1],  # G
    ], dtype=np.float32)

    fwd, rev = _scan_single_sequence(txt_seq, arr_motif)

    P = len(txt_seq) - arr_motif.shape[0] + 1
    assert fwd.shape == (P,)
    assert rev.shape == (P,)
    assert np.isfinite(fwd).all()
    assert np.isfinite(rev).all()


def test_scan_perfect_match_single_window():
    txt_seq = "ACGT"
    arr_motif = np.array([
        [2, 0, 0, 0],  # A
        [0, 3, 0, 0],  # C
        [0, 0, 5, 0],  # G
        [0, 0, 0, 7],  # T
    ], dtype=np.float32)

    fwd, rev = _scan_single_sequence(txt_seq, arr_motif)

    assert fwd.shape == (1,)
    assert np.isclose(fwd[0], 17.0)
    assert rev.shape == (1,)
    assert np.isfinite(rev[0])


def test_scan_single_mismatch():
    txt_seq = "ACGA"  # last base mismatched (should be T)
    arr_motif = np.array([
        [2, 0, 0, 0],  # A
        [0, 3, 0, 0],  # C
        [0, 0, 5, 0],  # G
        [0, 0, 0, 7],  # T
    ], dtype=np.float32)

    fwd, _ = _scan_single_sequence(txt_seq, arr_motif)

    assert fwd.shape == (1,)
    assert np.isclose(fwd[0], 10.0)  # 2+3+5+0


def test_scan_multiple_windows_values():
    txt_seq = "AACGT"
    arr_motif = np.array([
        [2, 0, 0, 0],  # A
        [0, 3, 0, 0],  # C
        [0, 0, 5, 0],  # G
    ], dtype=np.float32)

    fwd, _ = _scan_single_sequence(txt_seq, arr_motif)

    # windows: AAC, ACG, CGT
    expected = np.array([2+0+0, 2+3+5, 0+0+0], dtype=np.float32)
    assert np.allclose(fwd, expected)


def test_rc_consistency_forward_equals_rc_scanned_reversed():
    """
    score(seq, motif) == reverse(score(rc(seq), rc(motif)))
    """
    txt_seq = "ACGTACGT"
    arr_motif = np.array([
        [2, 0, 0, 0], # A
        [0, 3, 0, 0], # C
        [0, 0, 5, 0], # G
        [0, 0, 0, 7], # T
    ], dtype=np.float32)

    fwd1, _ = _scan_single_sequence(txt_seq, arr_motif)

    txt_seq_rc = reverse_complement_sequence(txt_seq)
    arr_motif_rc = reverse_complement_matrix(arr_motif, txt_alphabet="ACGT")

    fwd2, _ = _scan_single_sequence(txt_seq_rc, arr_motif_rc)

    assert np.allclose(fwd1, fwd2[::-1])


def test_sequence_too_short_raises():
    txt_seq = "ACG"
    arr_motif = np.ones((5, 4), dtype=np.float32)
    with pytest.raises(ValueError):
        _scan_single_sequence(txt_seq, arr_motif)
