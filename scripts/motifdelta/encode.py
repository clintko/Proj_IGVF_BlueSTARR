"""
Module for one-hot encoding sequence to matrix and decoding the corresponding matrix back to sequence
"""
import numpy as np

def one_hot_encode(txt_seq, txt_alphabet = "ACGT"):
    """
    one hot encoding of sequence (N,) to matrix (N,B) 
    where N is the sequence length and B is the alphabet size.

    Parameters
    ----------
    txt_seq : str
        DNA sequence (string of A/C/G/T).
    txt_alphabet : str
        order of base for one-hot encoding
        Default: "ACGT" (B=4)

    Returns
    -------
    np.ndarray
        matrix of one-hot encoded sequence, with shape (N,B)
        Default: shape = (N,4)
    """
    ### init array for sequence matrix
    arr = np.zeros((len(txt_seq), len(txt_alphabet)), dtype=np.float32)
    
    ### set base-position map for encoding sequence
    dct = {base: idx for idx, base in enumerate(txt_alphabet)}

    ### iterate through sequence and perform one-hot encoding
    for idx, base in enumerate(txt_seq):
        if base in dct:
            arr[idx, dct[base]] = 1.0
    return arr

def decode_one_hot(arr_seq_NxB, txt_alphabet = "ACGT"):
    """Helper function to turn one-hot encoded matrix back to seuqence"""
    return "".join(txt_alphabet[idx] for idx in arr_seq_NxB.argmax(axis=1))



