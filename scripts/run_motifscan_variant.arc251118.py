"""
Run motif scanning on sequences from FASTA file
"""

import numpy as np
import os, time
import argparse
import pickle

from Bio import SeqIO

#from fun_motifscan import scan_sequence_batch_parallel
from motifdelta import scan_sequence_batch_parallel


def main(args):
    # ==============================
    # Load FASTA sequences
    # ------------------------------
    lst_seq_record = list(SeqIO.parse(args.txt_fpath_fasta, "fasta"))
    print(f"Loaded {len(lst_seq_record)} sequences from {args.txt_fpath_fasta}\n")

    lst_txt_seq_idx = [rec.id for rec in lst_seq_record]
    lst_txt_seq_str = [str(rec.seq).upper() for rec in lst_seq_record]
    
    # ==============================
    # Load motif matrices
    # ------------------------------
    obj = np.load(args.txt_fpath_motif, allow_pickle=True)
    #dct_motif_lods = obj["lods"].item()
    #print(f"Loaded {len(dct_motif_lods)} motifs from {args.txt_fpath_motif}\n")

    dct_motif_model = obj["dct_results"].item()  # motif_name -> {arr_lod_WxB, grid, ccdf, Tbind, ...}
    print(f"Loaded PMAPs for {len(dct_motif_model)} motifs from {args.txt_fpath_motif}\n")
    
    # ==============================
    # Run motif scanning
    # ------------------------------
    print(f"Running motif scanning using {args.num_core} cores...")
    time_start = time.time()

    #dct_results = scan_sequence_batch_parallel(
    #    lst_seq,
    #    dct_motif_lods,
    #    num_workers=args.num_core
    #)
    dct_results = scan_sequence_batch_parallel(
        lst_txt_seq_idx,
        lst_txt_seq_str,
        dct_motif_model,
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
