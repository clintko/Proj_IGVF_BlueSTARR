"""
Run motif scanning on sequences from FASTA file
"""

import numpy as np
import pickle 
import time
import argparse

from Bio import SeqIO

from motifdelta.encode import one_hot_encode
from motifdelta.rc     import reverse_complement_matrix


def pad_motif_matrix(dct_motif_matrix):
    """
    Pad motif matrices to the same width (right-pad with zeros).
    
    Parameters
    ----------
    dct_motif_matrix : dict
        {motif_name: arr_matrix (W, B)}, e.g. LLR or PWM matrices.
    
    Returns
    -------
    arr_kernel : np.ndarray
        Stacked padded motif matrices, shape (M, W_max, B) float32,
        where M = number of motifs, W_max = max motif width.
    """
    ### get motif info
    lst_motif_name  = list(dct_motif_matrix.keys())
    lst_motif_width = [dct_motif_matrix[txt].shape[0] for txt in lst_motif_name]
    
    if len(lst_motif_name) == 0:
        raise ValueError("dct_motif_matrix is empty.")
        
    ### get motif parameters
    W_max = max(lst_motif_width)
    M     = len(lst_motif_name)
    B     = dct_motif_matrix[lst_motif_name[0]].shape[1]

    ### output motif kernels
    arr_motif_kernel = np.zeros((M, W_max, B), dtype=np.float32)
    for idx, txt in enumerate(lst_motif_name):
        
        ### get motif matrix
        arr = dct_motif_matrix[txt]
        W   = arr.shape[0]
        
        ### double check all motif matrix must have the same alphabet dimention
        assert arr.shape[1] == B, "All motif matrices must have the same B dimension"

        ### assign matrix to kernel
        arr_motif_kernel[idx, :W, :] = arr
        
    return arr_motif_kernel


def batch_sliding_window(arr_seq_NxLxB, W):
    """
    Create batched sliding windows over sequences.

    Parameters
    ----------
    arr_seq_NxLxB : np.ndarray
        Shape (N, L, B). One-hot encoded sequences (or any L×B features).
    W : int
        Window width (e.g. max motif width).

    Returns
    -------
    arr_out : np.ndarray
        Sliding windows of the one-hot encoded input sequence.
        Shape (N, P, W, B), where P = L - W + 1.
    """
    ### get dimensions
    N, L, B = arr_seq_NxLxB.shape
    P = L - W + 1
    
    if N <= 0:
        raise ValueError("arr_seq_NxLxB is empty.")
    if W <= 0:
        raise ValueError(f"Window width W must be positive (got W={W}).")
    if W > L:
        raise ValueError(f"Window width W={W} > sequence length L={L}.")

    ### initiate output: (N, P, W, B)
    arr_out = np.empty((N, P, W, B), dtype=arr_seq_NxLxB.dtype)

    ### Compute sliding windows per sequence
    for idx in range(N):
        ### create sliding windows; result shape: (P, 1, W, B)
        arr_win = np.lib.stride_tricks.sliding_window_view(
            arr_seq_NxLxB[idx], (W, B)
        )
        ### reshape and add to output
        arr_out[idx] = arr_win.reshape(P, W, B)

    return arr_out
    

def scan_sequence_batch(arr_seq_NxLxB, arr_motif_fwd_MxWxB, arr_motif_rev_MxWxB):
    """
    Batch motif scanning for forward and reverse strands.

    Parameters
    ----------
    arr_seq_NxLxB : np.ndarray
        One-hot encoded sequences, shape (N, L, B).
    arr_motif_fwd_MxWxB : np.ndarray
        Forward-strand motif kernels, shape (M, W, B).
    arr_motif_rev_MxWxB : np.ndarray
        Reverse-strand motif kernels, shape (M, W, B).

    Returns
    -------
    arr_out : np.ndarray
        Scanning results, shape (N, M, P, 2).
        2 = forward and reverse strands.
    """
    # ============================
    # Get dimensions
    # ----------------------------
    ### expect sequence and motif with three dimentions
    assert arr_seq_NxLxB.ndim       == 3, "Expected arr_seq_NxLxB with ndim=3"
    assert arr_motif_fwd_MxWxB.ndim == 3, "Expected arr_motif_fwd_MxWxB with ndim=3"
    assert arr_motif_rev_MxWxB.ndim == 3, "Expected arr_motif_rev_MxWxB with ndim=3"
    
    ### get dimensions
    N,  L,  B     = arr_seq_NxLxB.shape
    M,  W,  B_fwd = arr_motif_fwd_MxWxB.shape
    M2, W2, B_rev = arr_motif_rev_MxWxB.shape

    ### sanity check dimensions
    assert M == M2,  "Forward and reverse must have same number of motifs"
    assert W == W2,  "Forward and reverse must have same motif width"
    assert B == B_fwd == B_rev, "Alphabet dimension mismatch"

    ### calculate the dimention of score vector
    P = L - W + 1
    if P <= 0:
        raise ValueError(f"Sequence length L={L} must be >= motif width W={W}")

    # ============================
    # Scan motifs along sequence
    # ----------------------------
    
    ### create sliding windows
    ### arr_win: (N, P, W, B)
    arr_win = batch_sliding_window(arr_seq_NxLxB, W)

    ### forward strand scanning
    ### arr_fwd: (N, M, P)
    arr_fwd = np.einsum(
        "n p w b, m w b -> n m p",
        arr_win,
        arr_motif_fwd_MxWxB,
        optimize=True
    )

    ### reverse strand scanning
    ### arr_rev: (N, M, P)
    arr_rev = np.einsum(
        "n p w b, m w b -> n m p",
        arr_win,
        arr_motif_rev_MxWxB,
        optimize=True
    )

    ### combine forward + reverse
    ### arr_out: (N, M, P, 2)
    arr_out = np.stack([arr_fwd, arr_rev], axis=-1).astype(np.float32)
    return arr_out

    
def main(args):
    # ============================
    # Load FASTA sequences
    # ----------------------------

    print("Loading FASTA sequences...")
    time_start = time.time()
    
    ### import fasta as sequence record
    lst_seq_record_ref = list(SeqIO.parse(args.txt_fpath_fasta_ref, "fasta"))
    lst_seq_record_obs = list(SeqIO.parse(args.txt_fpath_fasta_obs, "fasta"))
    lst_seq_record_ubs = list(SeqIO.parse(args.txt_fpath_fasta_ubs, "fasta"))
    print(f"Loaded {len(lst_seq_record_ref)} sequences")

    if len(lst_seq_record_ref) == 0:
        raise ValueError("No sequences loaded from FASTA (ref).")

    ### convert Biopython records to strings
    lst_seq_txt_idx = [str(rec.id)          for rec in lst_seq_record_ref]
    lst_seq_txt_ref = [str(rec.seq).upper() for rec in lst_seq_record_ref]
    lst_seq_txt_obs = [str(rec.seq).upper() for rec in lst_seq_record_obs]
    lst_seq_txt_ubs = [str(rec.seq).upper() for rec in lst_seq_record_ubs]
    
    ### sanity check: sequence index
    assert [rec.id for rec in lst_seq_record_ref] == \
           [rec.id for rec in lst_seq_record_obs] == \
           [rec.id for rec in lst_seq_record_ubs], \
           "Ref/Obs/Ubs FASTA files do not share the same sequence IDs."

    ### sanity check: sequence length
    L = len(lst_seq_txt_ref[0])
    assert all(len(seq) == L for seq in lst_seq_txt_ref), "Ref sequences have inconsistent length."
    assert all(len(seq) == L for seq in lst_seq_txt_obs), "Obs sequences have inconsistent length."
    assert all(len(seq) == L for seq in lst_seq_txt_ubs), "Ubs sequences have inconsistent length."

    time_runtime = time.time() - time_start
    print(f"Load and check complete in {time_runtime:.2f} seconds\n")
    
    # ==============================
    # Load motif matrices
    # ------------------------------

    print("Loading motif matrices...")
    time_start = time.time()
    
    with open(args.txt_fpath_motif, "rb") as f:
        obj = pickle.load(f)

    dct_motif_lods = obj["lods"]
    txt_alphabet   = obj["alphabet"]
    print(f"Loaded {len(dct_motif_lods)} motifs from {args.txt_fpath_motif}")

    time_runtime = time.time() - time_start
    print(f"Load and check complete in {time_runtime:.2f} seconds\n")
    
    # ==============================
    # set motif kernel
    # ------------------------------

    print("Setting motif kernel...")
    time_start = time.time()
    
    ### build forward (unpadded -> right-padded)
    arr_motif_lods_fwd = pad_motif_matrix(dct_motif_lods)
    print("Forward kernels shape:", arr_motif_lods_fwd.shape)
    
    ### apply reverse complement on the unpadded motifs
    dct_motif_lods_rev = {
        txt: reverse_complement_matrix(arr) for txt, arr in dct_motif_lods.items()
    }
    
    ### right-pad the reverse-complement motifs in the same way
    arr_motif_lods_rev = pad_motif_matrix(dct_motif_lods_rev)
    print("Reverse kernels shape:", arr_motif_lods_rev.shape)

    #### get motif names
    lst_motif_name = list(dct_motif_lods.keys())

    time_runtime = time.time() - time_start
    print(f"Set complete in {time_runtime:.2f} seconds\n")
    
    # ==============================
    # Run motif scanning
    # ------------------------------
    print("Running motif scanning...")
    time_start = time.time()

    ### one-hot encoding sequences
    arr_seq_ref = np.stack([one_hot_encode(txt, txt_alphabet=txt_alphabet) for txt in lst_seq_txt_ref])
    arr_seq_obs = np.stack([one_hot_encode(txt, txt_alphabet=txt_alphabet) for txt in lst_seq_txt_obs])
    arr_seq_ubs = np.stack([one_hot_encode(txt, txt_alphabet=txt_alphabet) for txt in lst_seq_txt_ubs])

    ### motif scanning
    arr_motif_scan_ref = scan_sequence_batch(arr_seq_ref, arr_motif_lods_fwd, arr_motif_lods_rev)
    arr_motif_scan_obs = scan_sequence_batch(arr_seq_obs, arr_motif_lods_fwd, arr_motif_lods_rev)
    arr_motif_scan_ubs = scan_sequence_batch(arr_seq_ubs, arr_motif_lods_fwd, arr_motif_lods_rev)
    
    time_runtime = time.time() - time_start
    print(f"Scan complete in {time_runtime:.2f} seconds")
    
    num_memory = (arr_motif_scan_ref.nbytes * 3) / (1024**3)
    print(f"Estimated memory use (ref+obs+unobs): {num_memory:.3f} GB\n")
    
    # ==============================
    # Save results
    # ------------------------------
    
    print("Saving results...")
    time_start = time.time()
    
    #np.savez_compressed(
    np.savez(
        args.txt_fpath_output,
        Ref          = arr_motif_scan_ref,
        Obs          = arr_motif_scan_obs,
        Ubs          = arr_motif_scan_ubs,
        Idx_Sequence = np.array(lst_seq_txt_idx),
        Idx_Motif    = np.array(lst_motif_name)
    )

    time_runtime = time.time() - time_start
    
    print(f"Saved results to {args.txt_fpath_output}")
    print(f"Saved complete in {time_runtime:.2f} seconds\n")
    print("Done.")
    

if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files")
    
    parser.add_argument("--txt_fpath_fasta_ref",  type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--txt_fpath_fasta_obs",  type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--txt_fpath_fasta_ubs",  type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--txt_fpath_motif",      type=str, required=True, help="Path to motif pickle file")
    parser.add_argument("--txt_fpath_output",     type=str, required=True, help="Path to save motif scan results")

    args = parser.parse_args()

    ### run main function
    main(args)