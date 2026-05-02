"""
Run motif scanning on sequences from FASTA file
"""

import numpy as np
import pickle
import time
import argparse

from fun_fasta import fun_read_fasta

from motifdelta.score import (
    prepare_motif_kernels,
    one_hot_batch,
    scan_sequence_batch,
)


def parse_variant_id(txt_id: str):
    """
    Parse variant ID in the format:
        chr:pos:REF:OBS:UBS

    Returns
    -------
    tuple (chrom, pos, ref, obs, ubs)
    """
    parts = txt_id.split(":")
    if len(parts) != 5:
        raise ValueError(f"Unexpected variant id format (expected 5 fields): {txt_id}")
    chrom, pos, ref, obs, ubs = parts
    return chrom, int(pos), ref, obs, ubs


def replace_center_base(txt_seq: str, num_center: int, txt_alt: str) -> str:
    """
    Replace the base at num_center with txt_alt.

    Parameters
    ----------
    txt_seq : str
        Reference sequence window
    num_center : int
        0-based index of the variant position in the window
    txt_alt : str
        Alternate allele (single base)

    Returns
    -------
    str
        Mutated sequence window
    """
    if len(txt_alt) != 1:
        raise ValueError(f"Only SNVs supported here (got alt={txt_alt})")
    if not (0 <= num_center < len(txt_seq)):
        raise ValueError(f"Center index out of range: center={num_center}, L={len(txt_seq)}")
    return txt_seq[:num_center] + txt_alt + txt_seq[num_center + 1:]


def main(args):
    # ============================
    # Load FASTA sequences
    # ----------------------------
    
    print("Loading FASTA sequences...")
    time_start = time.time()

    ### import fasta as sequence record
    lst_seq_records = list(fun_read_fasta(args.txt_fpath_fasta))
    lst_txt_rec_idx = [idx for idx, _ in lst_seq_records]
    lst_txt_seq_ref = [seq.upper() for _, seq in lst_seq_records]
    print(f"Loaded {len(lst_seq_records)} sequences")

    if len(lst_seq_records) == 0:
        raise ValueError("No sequences loaded from FASTA (ref).")

    ### sanity check: sequence length
    num_length = len(lst_txt_seq_ref[0])
    if not all(len(seq) == num_length for seq in lst_txt_seq_ref):
        raise ValueError("Ref sequences have inconsistent length.")

    ### determine center index in the window
    ### For flankL=35, flankR=70, center = 35 (0-based)
    num_center = args.num_flank_left
    if (num_center < 0) or (num_center >= num_length):
        raise ValueError(f"Invalid num_flank_left={args.num_flank_left} for sequence length L={num_length}")

    ### generate Obs/Ubs sequences from Ref + variant id
    lst_txt_seq_obs = []
    lst_txt_seq_ubs = []

    for txt_rec_idx, txt_seq_ref in zip(lst_txt_rec_idx, lst_txt_seq_ref):
        ### parse variant id
        txt_chrom, num_pos, txt_allele_ref, txt_allele_obs, txt_allele_ubs = parse_variant_id(txt_rec_idx)

        ### sanity check: reference allele matches the ref sequence center (when A/C/G/T)
        if len(txt_allele_ref) != 1:
            raise ValueError(f"Only SNVs supported here (got ref={txt_allele_ref}) in {txt_rec_idx}")

        txt_base_center = txt_seq_ref[num_center]
        if (txt_base_center in "ACGT") and (txt_base_center != txt_allele_ref):
            raise ValueError(
                f"Ref allele mismatch at center for {txt_rec_idx}: "
                f"seq_center={txt_base_center}, header_ref={txt_allele_ref}"
            )

        ### generate Obs/Ubs sequences
        lst_txt_seq_obs.append(replace_center_base(txt_seq_ref, num_center, txt_allele_obs))
        lst_txt_seq_ubs.append(replace_center_base(txt_seq_ref, num_center, txt_allele_ubs))

    time_runtime = time.time() - time_start
    print(f"Load and check complete in {time_runtime:.2f} seconds\n")

    # ==============================
    # Load motif matrices
    # ------------------------------
    
    print("Loading motif matrices...")
    time_start = time.time()

    with open(args.txt_fpath_motif, "rb") as f:
        obj = pickle.load(f)

    ### get motif data/info
    dct_motif_lods = obj["lods"]
    txt_alphabet   = obj["alphabet"]
    lst_motif_name = obj["names"]

    ### sanity check: same length of motif name list
    if len(lst_motif_name) != len(dct_motif_lods):
        raise ValueError(f"names length {len(lst_motif_name)} != lods length {len(dct_motif_lods)}")
        
    ### sanity check: missing or extra motif in name list
    lst_motif_miss = [m for m in lst_motif_name if m not in dct_motif_lods]
    if lst_motif_miss:
        raise ValueError(
            f"{len(lst_motif_miss)} motif names in obj['names'] missing from obj['lods']; "
            f"example: {lst_motif_miss[:3]}"
        )
        
    set_motif_name = set(lst_motif_name)
    lst_motif_extra = [m for m in dct_motif_lods.keys() if m not in set_motif_name]
    if lst_motif_extra:
        raise ValueError(
            f"{len(lst_motif_extra)} motifs in obj['lods'] not present in obj['names']; "
            f"example: {lst_motif_extra[:3]}"
        )
        
    ### show progress
    print(f"Loaded {len(dct_motif_lods)} motifs from {args.txt_fpath_motif}")

    time_runtime = time.time() - time_start
    print(f"Load and check complete in {time_runtime:.2f} seconds\n")

    # ==============================
    # set motif kernel
    # ------------------------------
    
    print("Setting motif kernel...")
    time_start = time.time()

    ### build forward and reverse kernels (unpadded -> right-padded)
    arr_motif_lods_fwd, arr_motif_lods_rev, lst_name, W_max = prepare_motif_kernels(
        dct_motif_lods,
        lst_motif_name,
        txt_alphabet = txt_alphabet
    )
    print("Forward kernels shape:", arr_motif_lods_fwd.shape)
    print("Reverse kernels shape:", arr_motif_lods_rev.shape)
    if lst_name != lst_motif_name:
        raise RuntimeError("Motif name order changed unexpectedly.")
    
    time_runtime = time.time() - time_start
    print(f"Set complete in {time_runtime:.2f} seconds\n")
    
    # ==============================
    # Run motif scanning
    # ------------------------------
    
    print("Running motif scanning...")
    time_start = time.time()

    ### one-hot encoding sequences (Ref + generated Obs/Ubs)
    arr_seq_ref = one_hot_batch(lst_txt_seq_ref, txt_alphabet=txt_alphabet)
    arr_seq_obs = one_hot_batch(lst_txt_seq_obs, txt_alphabet=txt_alphabet)
    arr_seq_ubs = one_hot_batch(lst_txt_seq_ubs, txt_alphabet=txt_alphabet)

    ### motif scanning
    num_batch_size = None if args.batch_size == 0 else args.batch_size
    arr_motif_scan_ref = scan_sequence_batch(arr_seq_ref, arr_motif_lods_fwd, arr_motif_lods_rev, batch_size=num_batch_size)
    arr_motif_scan_obs = scan_sequence_batch(arr_seq_obs, arr_motif_lods_fwd, arr_motif_lods_rev, batch_size=num_batch_size)
    arr_motif_scan_ubs = scan_sequence_batch(arr_seq_ubs, arr_motif_lods_fwd, arr_motif_lods_rev, batch_size=num_batch_size)

    time_runtime = time.time() - time_start
    print(f"Scan complete in {time_runtime:.2f} seconds")

    num_memory = (arr_motif_scan_ref.nbytes * 3) / (1024**3)
    print(f"Output array size (ref+obs+unobs): {num_memory:.3f} GB\n")

    # ==============================
    # Save results
    # ------------------------------

    print("Saving results...")
    time_start = time.time()

    np.savez(
        args.txt_fpath_output,
        Ref = arr_motif_scan_ref.astype(np.float32, copy=False),
        Obs = arr_motif_scan_obs.astype(np.float32, copy=False),
        Ubs = arr_motif_scan_ubs.astype(np.float32, copy=False),
        Idx_Sequence = np.asarray(lst_txt_rec_idx, dtype=object),
        Idx_Motif    = np.asarray(lst_motif_name, dtype=object),
        Alphabet     = txt_alphabet
    )
    
    time_runtime = time.time() - time_start
    print(f"Saved results to {args.txt_fpath_output}")
    print(f"Saved complete in {time_runtime:.2f} seconds\n")

    
if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files")

    parser.add_argument("--txt_fpath_fasta",  type=str, required=True, help="Path to input FASTA (Ref) file")
    parser.add_argument("--txt_fpath_motif",  type=str, required=True, help="Path to motif pickle file")
    parser.add_argument("--txt_fpath_output", type=str, required=True, help="Path to save motif scan results")

    ### required to locate the variant in the window (0-based)
    parser.add_argument("--num_flank_left", type=int, default=35, help="Left flank length (center index). Default=35")
    parser.add_argument("--batch_size",     type=int, default=0,  help="Sequence batch size for scanning. 0 means scan all sequences at once.")
    
    args = parser.parse_args()

    ### run main function
    main(args)
