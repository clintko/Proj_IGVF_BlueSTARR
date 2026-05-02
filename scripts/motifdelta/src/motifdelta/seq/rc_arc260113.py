"""
Module for reverse compliment sequence and matrix
"""

import numpy as np


def reverse_complement_sequence(txt_seq):
    """
    Return the reverse complement of a DNA sequence string.

    Parameters
    ----------
    txt_seq : str
        DNA sequence (e.g. "ACGTN")

    Returns
    -------
    str
        Reverse complement of input sequence (e.g. "NACGT")
    """
    dct_complement = str.maketrans(
        "ACGTacgt",  ### from these alphabets
        "TGCAtgca"   ### to   those alphabets
    )
    return txt_seq.translate(dct_complement)[::-1] ### flip the order as well


def reverse_complement_matrix(arr_seq, txt_alphabet="ACGT"):
    """
    Reverse-complement one or more motif/sequence matrices.

    Automatically dispatches to single or batch mode depending
    on input dimensionality.

    Parameters
    ----------
    arr_seq : np.ndarray
        Array of shape (L,B) or (N,L,B)
    txt_alphabet : str
        Alphabet order, default "ACGT".

    Returns
    -------
    np.ndarray
        Reverse-complemented matrix/matrices of same shape.
    """
    arr = np.asarray(arr_seq)
    if arr.ndim == 2:
        return reverse_complement_matrix_single(arr, txt_alphabet=txt_alphabet)
    elif arr.ndim == 3:
        return reverse_complement_matrix_batch(arr, txt_alphabet=txt_alphabet)
    else:
        raise ValueError("Input must have shape (L,B) or (N,L,B)")


def reverse_complement_matrix_single(arr_seq_LxB, txt_alphabet = "ACGT"):
    """
    Reverse-complement a position × base matrix of shape (L,B) (Default: L,4).
    The input matrix reversed along positions and with columns 
    swapped A <-> T and C <-> G.

    Parameters
    ----------
    arr_seq_NxB : 
        one-hot enocded matrix of a sequence with shape (L, B).
        Default: shape (L, 4)

    Returns
    -------
    np.ndarray
        one-hot enocded matrix of the Reverse-complement sequence with shape (L, B)
        Default: shape (L, 4)
    """
    ### sanity check: alphabet must contain A, C, G, T
    required = set("ACGT")
    if not required.issubset(set(txt_alphabet)):
        raise ValueError(f"Alphabet must contain A,C,G,T. Got: {txt_alphabet}")
        
    ### map base -> column index
    dct_base_index = {base: idx for idx, base in enumerate(txt_alphabet)}
    
    ### define complement mapping
    dct_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    
    ### build index list for complement columns, preserving input order
    lst_swap_idx = [dct_base_index[dct_complement.get(base, base)] for base in txt_alphabet]

    ### flip positions + swap bases
    arr_seq_rc = arr_seq_LxB[::-1, lst_swap_idx]
    
    return arr_seq_rc

def reverse_complement_matrix_batch(arr_seq_NxLxB, txt_alphabet="ACGT"):
    """
    Reverse-complement one or more position × base matrices.

    Parameters
    ----------
    arr_seq_NxLxB : np.ndarray
        Array of shape (N, L, B) or (L, B), where
          N = number of sequence,
          L = sequence length,
          B = alphabet size (typically 4).
    txt_alphabet : str, optional
        Alphabet order (default "ACGT").
        Must contain at least A,C,G,T.

    Returns
    -------
    np.ndarray
        Reverse-complemented array with the same shape as input.

    Notes
    -----
    This function flips the position axis and swaps columns A<->T, C<->G.
    If input is 2D, it behaves identically to the reverse_complement_matrix().
    """
    required = set("ACGT")
    if not required.issubset(set(txt_alphabet)):
        raise ValueError(f"Alphabet must contain A,C,G,T. Got: {txt_alphabet}")

    ### map base -> column index
    dct_base_index = {base: idx for idx, base in enumerate(txt_alphabet)}
    
    ### define complement mapping
    dct_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    
    ### build index list for complement columns, preserving input order
    lst_swap_idx = [dct_base_index[dct_complement.get(base, base)] for base in txt_alphabet]

    ### flip positions + swap bases
    arr_seq_rc = arr_seq_NxLxB[:, ::-1, :][:, :, lst_swap_idx]
    
    return arr_seq_rc