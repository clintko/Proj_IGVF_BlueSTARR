"""
Module for reverse compliment sequence and matrix
"""

import numpy as np

_RC_IDX_CACHE: dict[str, np.ndarray] = {}


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


def _swap_idx_for_alphabet(txt_alphabet: str) -> np.ndarray:
    """
    Return column indices for complement swap in the given alphabet order.
    Cached by exact alphabet string.
    """
    idx = _RC_IDX_CACHE.get(txt_alphabet)
    if idx is not None:
        return idx

    # Validate alphabet contains A,C,G,T in any case
    if not set("ACGT").issubset(set(txt_alphabet.upper())):
        raise ValueError(f"Alphabet must contain A,C,G,T. Got: {txt_alphabet}")

    dct_base_index = {base: i for i, base in enumerate(txt_alphabet)}
    # complement mapping uses uppercase keys, but we preserve original-case bases in alphabet
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}

    swap = []
    for base in txt_alphabet:
        b_up = base.upper()
        b_comp = comp.get(b_up, b_up)
        # find a matching letter in the alphabet with same case as original base
        # if alphabet is uppercase, this is just b_comp
        # if alphabet is lowercase, use lowercase
        b_comp_same_case = b_comp.lower() if base.islower() else b_comp
        swap.append(dct_base_index[b_comp_same_case])

    idx = np.asarray(swap, dtype=np.int64)
    _RC_IDX_CACHE[txt_alphabet] = idx
    return idx


def reverse_complement_matrix(arr_seq: np.ndarray, txt_alphabet: str = "ACGT") -> np.ndarray:
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
    if arr.ndim == 3:
        return reverse_complement_matrix_batch(arr, txt_alphabet=txt_alphabet)
    raise ValueError("Input must have shape (L,B) or (N,L,B)")


def reverse_complement_matrix_single(arr_seq_LxB: np.ndarray, txt_alphabet: str = "ACGT") -> np.ndarray:
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
    swap_idx = _swap_idx_for_alphabet(txt_alphabet)
    return arr_seq_LxB[::-1, :][:, swap_idx]


def reverse_complement_matrix_batch(arr_seq_NxLxB: np.ndarray, txt_alphabet: str = "ACGT") -> np.ndarray:
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
    swap_idx = _swap_idx_for_alphabet(txt_alphabet)
    return arr_seq_NxLxB[:, ::-1, :][:, :, swap_idx]