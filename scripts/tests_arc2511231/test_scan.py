"""
Unit tests for motif scanning functions
"""

import numpy as np
import pytest

from motifdelta.rc   import reverse_complement_sequence, reverse_complement_matrix
from motifdelta.scan import scan_motif, scan_motif_both_strands


def test_scan_motif_basic():
    txt_seq = "ACGTAC"
    arr_motif = np.array([
        [ 1, -1, -1, -1],  # A
        [-1,  1, -1, -1],  # C
        [-1, -1,  1, -1]   # G
    ], dtype=float)

    scores = scan_motif(txt_seq, arr_motif)

    ### check finite values and shape for forward scan
    assert isinstance(scores, np.ndarray)
    assert scores.shape[0] == len(txt_seq) - arr_motif.shape[0] + 1
    assert np.isfinite(scores).all()


def test_scan_motif_reverse_complement():
    txt_seq = "ACGTAC"
    arr_motif = np.array([
        [ 1, -1, -1, -1],  # A
        [-1,  1, -1, -1],  # C
        [-1, -1,  1, -1]   # G
    ], dtype=float)

    arr_score_fwd, arr_score_rev = scan_motif_both_strands(txt_seq, arr_motif)

    ### check finite values and shape for forward/reverse scan
    assert isinstance(arr_score_fwd, np.ndarray)
    assert isinstance(arr_score_rev, np.ndarray)
    assert arr_score_fwd.shape == arr_score_rev.shape
    assert np.isfinite(arr_score_fwd).all()
    assert np.isfinite(arr_score_rev).all()


def test_scan_motif_perfect_match():
    txt_seq = "ACGT"
    arr_motif = np.array([
        [2, 0, 0, 0],  # A
        [0, 3, 0, 0],  # C
        [0, 0, 5, 0],  # G
        [0, 0, 0, 7],  # T
    ], dtype=float)

    arr_score = scan_motif(txt_seq, arr_motif)
    
    ### Expect one exact window, score = 2 + 3 + 5 + 7 = 17
    assert arr_score.shape == (1,)
    assert np.isclose(arr_score[0], 17.0)


def test_scan_motif_single_mismatch():
    txt_seq = "ACGA"  # last base mismatched (should be T)
    arr_motif = np.array([
        [2, 0, 0, 0],  # A
        [0, 3, 0, 0],  # C
        [0, 0, 5, 0],  # G
        [0, 0, 0, 7],  # T
    ], dtype=float)
    
    arr_score = scan_motif(txt_seq, arr_motif)
    
    ### Last base mismatch, expected score = 2 + 3 + 5 + 0 = 10
    assert arr_score.shape == (1,)
    assert np.isclose(arr_score[0], 10.0)


def test_scan_motif_multiple_windows():
    txt_seq = "AACGT"
    arr_motif = np.array([
        [2, 0, 0, 0],  # A
        [0, 3, 0, 0],  # C
        [0, 0, 5, 0],  # G
    ], dtype=float)
    
    arr_score = scan_motif(txt_seq, arr_motif)
    
    ### Sequence length 5, motif length 3 -> 5-3+1 = 3 windows
    ### windows: AAC, ACG, CGT
    ###     AAC: scores = 2+0+0 =  2
    ###     ACG: scores = 2+3+5 = 10
    ###     CGT: scores = 0+0+0 =  0
    expected = [2+0+0, 2+3+5, 0+0+0]
    assert np.allclose(arr_score, expected)


def test_scan_motif_reverse_complement_consistency():
    txt_seq = "ACGTACGT"
    arr_motif = np.array([
        [2, 0, 0, 0], # A
        [0, 3, 0, 0], # C
        [0, 0, 5, 0], # G
        [0, 0, 0, 7], # T
    ], dtype=float)

    txt_seq_rc   = reverse_complement_sequence(txt_seq)
    arr_motif_rc = reverse_complement_matrix(arr_motif)

    arr_score_fwd = scan_motif(txt_seq,    arr_motif)
    arr_score_rev = scan_motif(txt_seq_rc, arr_motif_rc)

    ### RC(seq) scanned with RC(motif) should match forward strand
    assert np.allclose(arr_score_fwd, arr_score_rev[::-1])


def test_scan_motif_sequence_too_short():
    txt_seq   = "ACG"
    arr_motif = np.ones((5, 4))
    with pytest.raises(ValueError):
        _ = scan_motif(txt_seq, arr_motif)
    