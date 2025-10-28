### set environment
import numpy as np
import os, time
import argparse

from Bio import SeqIO 

from fun_motifscan import (
    scan_sequence_batch_base,
    scan_sequence_batch_serial,
    scan_sequence_batch_parallel,
)


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
    # Benchmark batch sequence scan
    # ------------------------------

    txt_mode = args.txt_benchmark_mode.lower()
    num_seq  = args.num_test_seq
    num_core = args.num_core
    
    print(f"Running benchmark mode: {txt_mode}")
    print(f"Testing {num_seq} sequences with {num_core} workers\n")

    ### start the timer
    time_start = time.time()
    
    if txt_mode == "base":
        res = scan_sequence_batch_base(lst_seq[:num_seq], dct_motif_lods)
        
    elif txt_mode == "serial":
        res = scan_sequence_batch_serial(lst_seq[:num_seq], dct_motif_lods)
        
    elif txt_mode == "parallel":
        res = scan_sequence_batch_parallel(
            lst_seq[:num_seq], 
            dct_motif_lods, 
            num_workers=num_core
        )
        
    else:
        raise ValueError(f"Unknown benchmark mode: {txt_mode}")

    ### end the timer
    time_runtime = time.time() - time_start

    # ==============================
    # Show and save benchmark results
    # ------------------------------
    print(f"Benchmark complete:")
    print(f"Runtime ({txt_mode}): {time_runtime:.2f} seconds\n")
    
    ### save benchmark results
    with open(args.txt_fpath_output, "w") as f:
        f.write("Mode,Num_Seq,Num_Core,Time_Seconds\n")
        f.write(f"{txt_mode},{num_seq},{num_core},{time_runtime:.2f}\n")

    print(f"Results saved to {args.txt_fpath_output}\n")

if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Benchmark motif scanning performance")
    parser.add_argument("--txt_fpath_fasta",    type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--txt_fpath_motif",    type=str, required=True, help="Path to motif npz file")
    parser.add_argument("--txt_fpath_output",   type=str, required=True, help="Path to save benchmark runtime")
    parser.add_argument("--num_test_seq",       type=int, default=50,    help="Number of sequences to test")
    parser.add_argument("--num_core",           type=int, default=8,     help="Number of parallel workers")
    parser.add_argument("--txt_benchmark_mode", type=str, required=True, choices=["base", "serial", "parallel"], help="Benchmark mode: base, serial, or parallel")
    
    args = parser.parse_args()

    ### execute main function
    main(args)