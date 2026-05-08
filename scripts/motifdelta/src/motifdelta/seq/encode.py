"""
motifdelta.seq.encode

Fast one-hot encoding and decoding utilities.
"""

from __future__ import annotations
import numpy as np


_ASCII_MAP_CACHE: dict[str, np.ndarray] = {}

def _make_ascii_map(txt_alphabet: str) -> np.ndarray:
    table = _ASCII_MAP_CACHE.get(txt_alphabet)
    if table is not None:
        return table
    table = np.full(256, -1, dtype=np.int16)
    for i, ch in enumerate(txt_alphabet):
        table[ord(ch)] = i
        table[ord(ch.lower())] = i
    _ASCII_MAP_CACHE[txt_alphabet] = table
    return table


def one_hot_encode(txt_seq: str, txt_alphabet: str = "ACGT") -> np.ndarray:
    """
    One-hot encode a sequence into shape (L, B), dtype float32.
    Unknown characters (e.g. N) become all-zero rows.
    """
    L = len(txt_seq)
    B = len(txt_alphabet)
    out = np.zeros((L, B), dtype=np.float32)

    if L == 0:
        return out

    table = _make_ascii_map(txt_alphabet)

    ### Convert to uint8 bytes, map to indices
    s = txt_seq.encode("ascii", "replace")
    idx = table[np.frombuffer(s, dtype=np.uint8)]  # (L,)

    valid = idx >= 0
    if np.any(valid):
        rows = np.nonzero(valid)[0]
        cols = idx[valid].astype(np.int64, copy=False)
        out[rows, cols] = 1.0

    return out


def decode_one_hot(arr_seq_LxB: np.ndarray, txt_alphabet: str = "ACGT") -> str:
    """Turn one-hot encoded matrix back to sequence using argmax."""
    if arr_seq_LxB.size == 0:
        return ""
    return "".join(txt_alphabet[i] for i in arr_seq_LxB.argmax(axis=1))


def one_hot_encode_batch(lst_seq: list[str], txt_alphabet: str = "ACGT") -> np.ndarray:
    """
    Encode a batch of equal-length sequences into (N, L, B).
    Unknown chars map to all-zero rows.
    """
    if len(lst_seq) == 0:
        raise ValueError("lst_seq is empty.")
    L0 = len(lst_seq[0])
    if not all(len(s) == L0 for s in lst_seq):
        raise ValueError("All sequences must have the same length.")

    table = _make_ascii_map(txt_alphabet)
    N = len(lst_seq)
    B = len(txt_alphabet)

    out = np.zeros((N, L0, B), dtype=np.float32)

    ### build a (N, L) uint8 array of bases
    buf = np.frombuffer(("".join(lst_seq)).encode("ascii", "replace"), dtype=np.uint8).reshape(N, L0)
    idx = table[buf]  # (N, L)

    valid = idx >= 0
    if np.any(valid):
        n_idx, l_idx = np.nonzero(valid)
        b_idx = idx[valid].astype(np.int64, copy=False)
        out[n_idx, l_idx, b_idx] = 1.0

    return out
