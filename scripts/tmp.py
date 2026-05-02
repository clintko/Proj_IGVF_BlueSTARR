"""
Run motif scanning on sequences from FASTA file
"""

import numpy as np
import pickle
import time
import argparse

from Bio import SeqIO

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
    Replace the base at num_center with txt_alt (SNV-only).

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
    lst_seq_record_ref = list(SeqIO.parse(args.txt_fpath_fasta_ref, "fasta"))
    print(f"Loaded {len(lst_seq_record_ref)} sequences")

    if len(lst_seq_record_ref) == 0:
        raise ValueError("No sequences loaded from FASTA (ref).")

    ### convert Biopython records to strings
    lst_seq_txt_idx = [str(rec.id)          for rec in lst_seq_record_ref]
    lst_seq_txt_ref = [str(rec.seq).upper() for rec in lst_seq_record_ref]

    ### sanity check: sequence length
    L = len(lst_seq_txt_ref[0])
    assert all(len(seq) == L for seq in lst_seq_txt_ref), "Ref sequences have inconsistent length."

    ### determine center index in the window
    ### For flankL=35, flankR=70, center = 35 (0-based)
    num_center = args.num_flank_left
    if (num_center < 0) or (num_center >= L):
        raise ValueError(f"Invalid num_flank_left={args.num_flank_left} for sequence length L={L}")

    ### generate Obs/Ubs sequences from Ref + variant id
    lst_seq_txt_obs = []
    lst_seq_txt_ubs = []

    for txt_id, txt_ref_seq in zip(lst_seq_txt_idx, lst_seq_txt_ref):
        _, _, txt_ref, txt_obs, txt_ubs = parse_variant_id(txt_id)

        ### sanity check: reference allele matches the ref sequence center (when A/C/G/T)
        if len(txt_ref) != 1:
            raise ValueError(f"Only SNVs supported here (got ref={txt_ref})")

        base_center = txt_ref_seq[num_center]
        if base_center in "ACGT" and base_center != txt_ref:
            raise ValueError(
                f"Ref allele mismatch at center for {txt_id}: "
                f"seq_center={base_center}, header_ref={txt_ref}"
            )

        ### generate
        lst_seq_txt_obs.append(replace_center_base(txt_ref_seq, num_center, txt_obs))
        lst_seq_txt_ubs.append(replace_center_base(txt_ref_seq, num_center, txt_ubs))

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

    ### build forward and reverse kernels (unpadded -> right-padded)
    arr_motif_lods_fwd, arr_motif_lods_rev, lst_motif_name, W_max = prepare_motif_kernels(
        dct_motif_lods,
        txt_alphabet=txt_alphabet
    )
    print("Forward kernels shape:", arr_motif_lods_fwd.shape)
    print("Reverse kernels shape:", arr_motif_lods_rev.shape)

    time_runtime = time.time() - time_start
    print(f"Set complete in {time_runtime:.2f} seconds\n")

    # ==============================
    # Run motif scanning
    # ------------------------------
    print("Running motif scanning...")
    time_start = time.time()

    ### one-hot encoding sequences (Ref + generated Obs/Ubs)
    arr_seq_ref = one_hot_batch(lst_seq_txt_ref, txt_alphabet=txt_alphabet)
    arr_seq_obs = one_hot_batch(lst_seq_txt_obs, txt_alphabet=txt_alphabet)
    arr_seq_ubs = one_hot_batch(lst_seq_txt_ubs, txt_alphabet=txt_alphabet)

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
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files (Ref-only FASTA)")

    parser.add_argument("--txt_fpath_fasta_ref", type=str, required=True, help="Path to input FASTA file (Ref only)")
    parser.add_argument("--txt_fpath_motif",     type=str, required=True, help="Path to motif pickle file")
    parser.add_argument("--txt_fpath_output",    type=str, required=True, help="Path to save motif scan results")

    ### required to locate the variant in the window (0-based)
    parser.add_argument("--num_flank_left", type=int, default=35, help="Left flank length (center index). Default=35")

    args = parser.parse_args()

    ### run main function
    main(args)
