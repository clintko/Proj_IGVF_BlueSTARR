"""
Module for estimating background nucleotide models
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