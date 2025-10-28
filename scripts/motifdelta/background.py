"""
Module for estimating background nucleotide models (0th- and 1st-order)
"""

import numpy as np
from collections import Counter


def background_zero_order(lst_txt_seq, alphabet="ACGT", pseudocount=1e-6):
    """
    Estimate 0th-order (i.i.d.) background nucleotide frequencies.

    Parameters
    ----------
    lst_txt_seq : list[str]
        List of DNA sequences (strings of ACGT)
    alphabet : str
        Order of bases (default: ACGT)
    pseudocount : float
        Small value added to avoid zeros

    Returns
    -------
    np.ndarray
        Array of background probabilities of length len(alphabet)
    """
    counts = Counter()
    total = 0

    for txt_seq in lst_txt_seq:
        for base in txt_seq:
            if base in alphabet:
                counts[base] += 1
                total += 1

    arr = np.array([counts.get(base, 0) for base in alphabet], dtype=float)
    arr = arr + pseudocount
    arr /= arr.sum()
    return arr


def background_first_order(lst_txt_seq, alphabet="ACGT", pseudocount=1e-6):
    """
    Estimate 1st-order Markov background: transition matrix P(x -> y).

    Parameters
    ----------
    lst_txt_seq : list[str]
        List of DNA sequences (strings of ACGT)
    alphabet : str
        Order of bases (default: ACGT)
    pseudocount : float
        Small value added to avoid zeros

    Returns
    -------
    np.ndarray
        Matrix shape (B,B) where entry [i,j] = P(y=j | x=i); B := len(alphabet)
    """
    B = len(alphabet)
    counts = np.zeros((B, B), dtype=float)
    dct_idx = {b: i for i, b in enumerate(alphabet)}

    for txt_seq in lst_txt_seq:
        for i in range(len(txt_seq) - 1):
            a, b = txt_seq[i], txt_seq[i + 1]
            if a in dct_idx and b in dct_idx:
                counts[dct_idx[a], dct_idx[b]] += 1

    counts += pseudocount
    row_sums = counts.sum(axis=1, keepdims=True)
    np.divide(counts, row_sums, out=counts, where=row_sums != 0)
    return counts