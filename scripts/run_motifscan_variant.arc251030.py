"""
Run motif scanning on sequences from FASTA file
"""

import numpy as np
import os, time
import argparse
import pickle

from Bio import SeqIO

from fun_motifscan import scan_sequence_batch_parallel


def main(args):
    # ==============================
    # Load FASTA sequences
    # ------------------------------
    lst_seq = list(SeqIO.parse(args.txt_fpath_fasta, "fasta"))
    print(f"Loaded {len(lst_seq)} sequences from {args.txt_fpath_fasta}\n")
    
    # ==============================
    # Load motif matrices
    # ------------------------------
    obj = np.load(args.txt_fpath_motif, allow_pickle=True)
    dct_motif_lods = obj["lods"].item()
    print(f"Loaded {len(dct_motif_lods)} motifs from {args.txt_fpath_motif}\n")

    # ==============================
    # Run motif scanning
    # ------------------------------
    print(f"Running motif scanning using {args.num_core} cores...")
    time_start = time.time()

    dct_results = scan_sequence_batch_parallel(
        lst_seq,
        dct_motif_lods,
        num_workers=args.num_core
    )

    time_runtime = time.time() - time_start
    print(f"Scan complete in {time_runtime:.2f} seconds\n")

    # ==============================
    # Save results
    # ------------------------------
    with open(args.txt_fpath_output, "wb") as file:
        pickle.dump(dct_results, file)

    print(f"Saved results to {args.txt_fpath_output}")


if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files")
    parser.add_argument("--txt_fpath_fasta",  type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--txt_fpath_motif",  type=str, required=True, help="Path to motif npz file")
    parser.add_argument("--txt_fpath_output", type=str, required=True, help="Path to save motif scan results")
    parser.add_argument("--num_core",         type=int, default=10,    help="Number of parallel cores")

    args = parser.parse_args()

    ### run main function
    main(args)
