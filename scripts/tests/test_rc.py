"""
Unit tests for reverse complement functions
"""

import numpy as np
from motifdelta.encode import one_hot_encode, decode_one_hot
from motifdelta.rc import reverse_complement_sequence, reverse_complement_matrix


def test_reverse_complement_sequence_basic():
    txt_seq_inp = "ACGTT"
    txt_seq_ans = "AACGT"
    txt_seq_out = reverse_complement_sequence(txt_seq_inp)
    assert txt_seq_out == txt_seq_ans


def test_reverse_complement_sequence_case_insensitive():
    txt_seq_inp = "aCgT"
    txt_seq_ans = "AcGt"
    txt_seq_out = reverse_complement_sequence(txt_seq_inp)
    assert txt_seq_out == txt_seq_ans


def test_reverse_complement_matrix_matches_string():
    txt_seq_inp = "ACGTT"
    txt_seq_ans = "AACGT"

    arr_seq_inp = one_hot_encode(txt_seq_inp)
    arr_seq_out = reverse_complement_matrix(arr_seq_inp, "ACGT")
    txt_seq_out = decode_one_hot(arr_seq_out)

    assert txt_seq_out == txt_seq_ans


def test_reverse_complement_matrix_shape():
    txt_seq_inp = "ACGTT"
    arr_seq_inp = one_hot_encode(txt_seq_inp)
    arr_seq_out = reverse_complement_matrix(arr_seq_inp, "ACGT")

    assert arr_seq_inp.shape == arr_seq_out.shape


def test_reverse_complement_matrix_double_rc_is_identity():
    txt_seq_inp = "ATGCATGC"
    arr_seq_inp = one_hot_encode(txt_seq_inp)
    arr_seq_tmp = reverse_complement_matrix(arr_seq_inp)
    arr_seq_out = reverse_complement_matrix(arr_seq_tmp)

    # Double reverse complement should return original
    assert np.allclose(arr_seq_inp, arr_seq_out)
    assert decode_one_hot(arr_seq_inp) == decode_one_hot(arr_seq_out)
