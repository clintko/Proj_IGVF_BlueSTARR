"""
Unit tests for one-hot encoding and decoding functions
"""

import numpy as np
from motifdelta.encode import one_hot_encode, decode_one_hot

def test_one_hot_case01_encode():
    txt_inp = "ACGTT"
    arr_inp = one_hot_encode(txt_inp)
    
    # Expected: each row corresponds to A,C,G,T,T
    arr_ans = np.array([
        [1,0,0,0],  # A
        [0,1,0,0],  # C
        [0,0,1,0],  # G
        [0,0,0,1],  # T
        [0,0,0,1],  # T
    ], dtype=np.float32)

    assert arr_inp.shape == (5,4)
    assert np.array_equal(arr_inp, arr_ans)


def test_one_hot_case02_decode():
    txt_inp = "ACGTT"
    arr_inp = one_hot_encode(txt_inp)
    txt_out = decode_one_hot(arr_inp)
    
    assert arr_inp.shape == (5,4)
    assert txt_inp == txt_out


def test_one_hot_case03_unknown_bases():
    txt_seq = "ANCT"

    # "N" should be zeros if alphabet = "ACGT"
    txt_alphabet = "ACGT"
    arr_seq = one_hot_encode(txt_seq, txt_alphabet=txt_alphabet)
    assert (arr_seq[1] == 0).all()

    # "N" should be [0,0,0,0,1] if alphabet = "ACGTN"
    txt_alphabet = "ACGTN"
    arr_seq = one_hot_encode(txt_seq, txt_alphabet=txt_alphabet)
    assert np.array_equal(arr_seq[1], np.array([0,0,0,0,1], dtype=np.float32))

    # "N" should be [0,1,0,0,0] if alphabet = "ANCGT"
    txt_alphabet = "ANCGT"
    arr_seq = one_hot_encode(txt_seq, txt_alphabet=txt_alphabet)
    assert np.array_equal(arr_seq[1], np.array([0,1,0,0,0], dtype=np.float32))