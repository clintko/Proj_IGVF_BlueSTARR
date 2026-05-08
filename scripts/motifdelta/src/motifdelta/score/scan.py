"""
motifdelta.score.scan

Fast batched motif scanning using padded motif kernels + einsum.

Design:
- Library code is array-only (no FASTA / no file I/O).
- Pipeline scripts handle reading FASTA and saving outputs.
"""

from __future__ import annotations

from typing import List, Tuple, Mapping

import numpy as np

from motifdelta.seq.encode import one_hot_encode
from motifdelta.seq.rc import reverse_complement_matrix


def pad_motif_matrix(
    dct_motif_matrix: Mapping[str, np.ndarray],
    lst_motif_name:   List[str] | None = None
) -> Tuple[np.ndarray, List[str], int]:
    """
    Pad motif matrices to the same width (right-pad with zeros).

    Parameters
    ----------
    dct_motif_matrix : dict
        {motif_name: arr_matrix (W, B)}.
    lst_motif_name : list[str]
        Motif names in the same order as kernel axis 0
        
    Returns
    -------
    arr_motif_kernel : np.ndarray
        Stacked padded motif matrices, shape (M, W_max, B) float32
    lst_motif_name : list[str]
        Motif name order in the matrix
    W_max : int
        Maximum motif width
    """
    ### set lst_motif_name variable
    if lst_motif_name is None:
        lst_motif_name = list(dct_motif_matrix.keys())

    ### validate
    if len(lst_motif_name) == 0:
        raise ValueError("No motifs.")
    if len(lst_motif_name) != len(set(lst_motif_name)):
        raise ValueError("Duplicate motif names in lst_motif_name.")
    missing = [m for m in lst_motif_name if m not in dct_motif_matrix]
    if missing:
        raise ValueError(f"Motifs missing from dictionary: {missing[:5]} (n={len(missing)})")

    ### get motif dimention
    lst_motif_width = [dct_motif_matrix[m].shape[0] for m in lst_motif_name]
    W_max = int(max(lst_motif_width))
    M = len(lst_motif_name)
    B = int(dct_motif_matrix[lst_motif_name[0]].shape[1])

    ### init: motif kernel
    arr_motif_kernel = np.zeros((M, W_max, B), dtype=np.float32)

    ### loop through motifs and set kernel for each motif
    for i, m in enumerate(lst_motif_name):
        ### get motif kernel
        arr = np.asarray(dct_motif_matrix[m], dtype=np.float32)

        ### validate 
        if arr.ndim != 2:
            raise ValueError(f"Motif {m} matrix must be 2D (W,B); got shape {arr.shape}")
        if arr.shape[1] != B:
            raise ValueError(f"Motif {m} has B={arr.shape[1]} but expected B={B}")

        ### set motif kernel
        W = int(arr.shape[0])
        arr_motif_kernel[i, :W, :] = arr

    return arr_motif_kernel, lst_motif_name, W_max
    

def prepare_motif_kernels(
    dct_motif_WxB: Mapping[str, np.ndarray],
    lst_motif_name: List[str] | None = None,
    txt_alphabet: str = "ACGT",
) -> Tuple[np.ndarray, np.ndarray, List[str], int]:
    """
    Build padded forward and reverse-complement motif kernels.

    Parameters
    ----------
    dct_motif_WxB : Mapping[str, np.ndarray]
        {motif_name: arr_matrix (W, B)}.
    lst_motif_name : list[str] | None
        Motif names in the desired order (kernel axis 0). If None, uses dict key order.
    txt_alphabet : str
        Alphabet for reverse-complement.

    Returns
    -------
    arr_fwd_MxWxB : np.ndarray
    arr_rev_MxWxB : np.ndarray
    lst_motif_name : list[str]
    W_max : int
    """
    ### forward kernel in canonical order
    arr_fwd, lst_name, W_max = pad_motif_matrix(dct_motif_WxB, lst_motif_name)

    ### reverse kernels computed in the SAME canonical order
    dct_rev = {
        m: reverse_complement_matrix(
            np.asarray(dct_motif_WxB[m], dtype=np.float32),
            txt_alphabet=txt_alphabet
        )
        for m in lst_name
    }
    arr_rev, lst_name2, W_max2 = pad_motif_matrix(dct_rev, lst_name)

    ### validate
    if lst_name2 != lst_name:
        raise RuntimeError("Motif name order mismatch between forward and reverse kernels.")
    if W_max2 != W_max:
        raise RuntimeError("W_max mismatch between forward and reverse kernels.")

    return arr_fwd, arr_rev, lst_name, W_max


def batch_sliding_window(arr_seq_NxLxB: np.ndarray, W: int) -> np.ndarray:
    """
    Create batched sliding windows over one-hot sequences.

    Input
    -----
    arr_seq_NxLxB : (N, L, B)

    Output
    ------
    arr_win : (N, P, W, B) where P = L - W + 1
    """
    if arr_seq_NxLxB.ndim != 3:
        raise ValueError(f"Expected (N,L,B) array; got shape {arr_seq_NxLxB.shape}")

    N, L, B = arr_seq_NxLxB.shape
    if N <= 0:
        raise ValueError("arr_seq_NxLxB is empty.")
    if W <= 0:
        raise ValueError(f"W must be positive (got W={W}).")
    if W > L:
        raise ValueError(f"W={W} > sequence length L={L}.")

    P = L - W + 1
    arr_out = np.empty((N, P, W, B), dtype=arr_seq_NxLxB.dtype)

    # per-sequence sliding_window_view (keeps memory predictable; avoids tricky strides across N)
    for i in range(N):
        win = np.lib.stride_tricks.sliding_window_view(arr_seq_NxLxB[i], (W, B))
        # win is (P, 1, W, B)
        arr_out[i] = win.reshape(P, W, B)

    return arr_out


def scan_sequence_batch(
    arr_seq_NxLxB: np.ndarray,
    arr_motif_fwd_MxWxB: np.ndarray,
    arr_motif_rev_MxWxB: np.ndarray,
    batch_size: int | None = None,
) -> np.ndarray:
    """
    Batch motif scanning for forward and reverse strands.

    Parameters
    ----------
    arr_seq_NxLxB : np.ndarray
        One-hot sequences, shape (N, L, B)
    arr_motif_fwd_MxWxB : np.ndarray
        Forward motif kernels, shape (M, W, B)
    arr_motif_rev_MxWxB : np.ndarray
        Reverse motif kernels, shape (M, W, B)
    batch_size : int | None
        If provided, scan sequences in chunks of this size to limit peak memory.

    Returns
    -------
    arr_out : np.ndarray
        Shape (N, M, P, 2) with forward/reverse in last dim.
    """
    ### sanity check data dimensions
    if arr_seq_NxLxB.ndim != 3:
        raise ValueError("Expected arr_seq_NxLxB with ndim=3")
    if arr_motif_fwd_MxWxB.ndim != 3:
        raise ValueError("Expected arr_motif_fwd_MxWxB with ndim=3")
    if arr_motif_rev_MxWxB.ndim != 3:
        raise ValueError("Expected arr_motif_rev_MxWxB with ndim=3")

    ### get data dimensions
    N, L, B = arr_seq_NxLxB.shape
    M, W, B1 = arr_motif_fwd_MxWxB.shape
    M2, W2, B2 = arr_motif_rev_MxWxB.shape

    if (M != M2) or (W != W2) or (B != B1) or (B != B2):
        raise ValueError(
            f"Dimension mismatch: seq (N,L,B)=({N},{L},{B}), "
            f"motif_fwd (M,W,B)=({M},{W},{B1}), motif_rev (M,W,B)=({M2},{W2},{B2})"
        )

    P = L - W + 1
    if P <= 0:
        raise ValueError(f"Sequence length L={L} must be >= motif width W={W}")


    ### decide batch size
    if batch_size is None:
        batch_size = N
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive (got {batch_size})")

    ### preallocate outputs (to avoid repeated large allocations)
    out_fwd = np.empty((N, M, P), dtype=np.float32)
    out_rev = np.empty((N, M, P), dtype=np.float32)

    ### loop over batches of sequences
    for i0 in range(0, N, batch_size):
        i1 = min(i0 + batch_size, N)
        arr_seq_batch = arr_seq_NxLxB[i0:i1]  # (Nb, L, B)

        # (Nb, P, W, B)
        arr_win = batch_sliding_window(arr_seq_batch, W)

        # (Nb, M, P)
        arr_fwd = np.einsum("n p w b, m w b -> n m p", arr_win, arr_motif_fwd_MxWxB, optimize=True)
        arr_rev = np.einsum("n p w b, m w b -> n m p", arr_win, arr_motif_rev_MxWxB, optimize=True)

        out_fwd[i0:i1] = arr_fwd.astype(np.float32, copy=False)
        out_rev[i0:i1] = arr_rev.astype(np.float32, copy=False)

    return np.stack([out_fwd, out_rev], axis=-1)


def one_hot_batch(lst_seq: List[str], txt_alphabet: str = "ACGT") -> np.ndarray:
    """
    One-hot encode a batch of sequences into (N,L,B).
    Assumes all sequences same length.
    """
    if len(lst_seq) == 0:
        raise ValueError("lst_seq is empty.")
    L0 = len(lst_seq[0])
    if not all(len(s) == L0 for s in lst_seq):
        raise ValueError("All sequences must have the same length for batch encoding.")
    return np.stack([one_hot_encode(s, txt_alphabet=txt_alphabet) for s in lst_seq])
