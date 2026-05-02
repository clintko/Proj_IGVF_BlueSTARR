"""
Compute motif delta arrays and motif gain/loss based on Tbind.

Input:
    <chunk>_scan.npz
    motif_tbind.pkl  (dict: {"meta": ..., "tbind": {motif_name: num_Tbind}})

Output (reduced):
    <chunk>_delta_reduced.npz
    <chunk>_summary_motif.tsv
    <chunk>_summary_variant.tsv
    <chunk>_event.tsv  (combined gain + loss events)

Notes
-----
This version avoids saving the full (N, M, P, 2) Delta/Gain/Loss arrays.
Instead, it reduces to at most one Gain and one Loss record per (variant, motif).
"""

import numpy as np
import pandas as pd
import argparse
import pickle
import time
import os


def compute_gain_loss(arr_obs, arr_ubs, arr_Tbind):
    """
    Compute boolean gain/loss masks based on Tbind thresholds.

    Parameters
    ----------
    arr_obs : np.ndarray
        Observed-allele scores, shape (N, M, P, 2).
    arr_ubs : np.ndarray
        Unobserved-allele scores, shape (N, M, P, 2).
    arr_Tbind : np.ndarray
        Motif thresholds, shape (M,).

    Returns
    -------
    arr_gain : np.ndarray (bool)
        True where motif is formed in Ubs but not in Obs, shape (N, M, P, 2).
    arr_loss : np.ndarray (bool)
        True where motif is present in Obs but not in Ubs, shape (N, M, P, 2).
    """
    ### dimention check
    if arr_obs.shape != arr_ubs.shape:
        raise ValueError(f"Shape mismatch between Obs {arr_obs.shape} and Ubs {arr_ubs.shape}")

    if arr_Tbind.ndim != 1:
        raise ValueError(f"arr_Tbind must be 1D, got ndim={arr_Tbind.ndim}")

    ### get dimention
    N, M, P, S = arr_obs.shape
    if arr_Tbind.shape[0] != M:
        raise ValueError(
            f"Tbind length M={arr_Tbind.shape[0]} does not match motif dimension M={M}"
        )

    ### broadcasting: (M,) -> (1, M, 1, 1)
    arr_Tbind_broadcast = arr_Tbind[None, :, None, None]

    ### calculate masks
    arr_gain = (arr_ubs >= arr_Tbind_broadcast) & (arr_obs <  arr_Tbind_broadcast)
    arr_loss = (arr_ubs <  arr_Tbind_broadcast) & (arr_obs >= arr_Tbind_broadcast)

    return arr_gain, arr_loss


def reduce_best_gain_loss(arr_obs, arr_ubs, arr_Tbind):
    """
    Reduce (N, M, P, 2) arrays into per-(variant,motif) best gain and best loss.

    We keep:
    - max gain delta (and where it occurs)
    - min loss delta (and where it occurs)

    Returns
    -------
    dct_out : dict[str, np.ndarray]
        Keys include:
            Gain_Any   (N,M) bool
            Loss_Any   (N,M) bool
            Gain_Delta (N,M) float32
            Loss_Delta (N,M) float32
            Gain_Index (N,M) int32   index in flattened (P*2)
            Loss_Index (N,M) int32
            Gain_Obs   (N,M) float32  Obs score at best gain site
            Gain_Ubs   (N,M) float32  Ubs score at best gain site
            Loss_Obs   (N,M) float32  Obs score at best loss site
            Loss_Ubs   (N,M) float32  Ubs score at best loss site
    """
    # ================================
    # Compute delta + gain/loss masks
    # --------------------------------

    ### delta
    arr_delta = arr_ubs - arr_obs  # (N,M,P,2)

    ### gain/loss masks
    arr_gain, arr_loss = compute_gain_loss(arr_obs, arr_ubs, arr_Tbind)

    # ================================
    # Reduce over (P, strand)
    # --------------------------------
    N, M, P, S = arr_obs.shape
    if S != 2:
        raise ValueError(f"Expected 2 strands, got S={S}")

    ### flatten last dims: (N,M,P,2) -> (N,M,P*2)
    arr_delta_flat = arr_delta.reshape(N, M, -1)
    arr_obs_flat   = arr_obs.reshape(N, M, -1)
    arr_ubs_flat   = arr_ubs.reshape(N, M, -1)
    arr_gain_flat  = arr_gain.reshape(N, M, -1)
    arr_loss_flat  = arr_loss.reshape(N, M, -1)

    ### any gain/loss per (N,M)
    arr_gain_any = arr_gain_flat.any(axis=2)  # (N,M)
    arr_loss_any = arr_loss_flat.any(axis=2)  # (N,M)

    # --------------------------------
    # Best gain: max delta among gain sites
    # --------------------------------
    ### -inf where not gain so max ignores them
    arr_delta_gain = np.where(arr_gain_flat, arr_delta_flat, -np.inf)  # (N,M,P*2)
    arr_gain_idx   = arr_delta_gain.argmax(axis=2).astype(np.int32)    # (N,M)
    arr_gain_delta = arr_delta_gain.max(axis=2).astype(np.float32)     # (N,M)

    ### gather obs/ubs at best gain site
    arr_gain_obs = np.take_along_axis(arr_obs_flat, arr_gain_idx[:, :, None], axis=2)[:, :, 0].astype(np.float32)
    arr_gain_ubs = np.take_along_axis(arr_ubs_flat, arr_gain_idx[:, :, None], axis=2)[:, :, 0].astype(np.float32)

    ### clean invalid cells (no gain): set delta to 0 and idx to -1, obs/ubs to nan
    arr_gain_delta = np.where(arr_gain_any, arr_gain_delta, 0.0).astype(np.float32)
    arr_gain_idx   = np.where(arr_gain_any, arr_gain_idx, -1).astype(np.int32)
    arr_gain_obs   = np.where(arr_gain_any, arr_gain_obs, np.nan).astype(np.float32)
    arr_gain_ubs   = np.where(arr_gain_any, arr_gain_ubs, np.nan).astype(np.float32)

    # --------------------------------
    # Best loss: min delta among loss sites
    # --------------------------------
    ### +inf where not loss so min ignores them
    arr_delta_loss = np.where(arr_loss_flat, arr_delta_flat, np.inf)   # (N,M,P*2)
    arr_loss_idx   = arr_delta_loss.argmin(axis=2).astype(np.int32)    # (N,M)
    arr_loss_delta = arr_delta_loss.min(axis=2).astype(np.float32)     # (N,M)

    ### gather obs/ubs at best loss site
    arr_loss_obs = np.take_along_axis(arr_obs_flat, arr_loss_idx[:, :, None], axis=2)[:, :, 0].astype(np.float32)
    arr_loss_ubs = np.take_along_axis(arr_ubs_flat, arr_loss_idx[:, :, None], axis=2)[:, :, 0].astype(np.float32)

    ### clean invalid cells (no loss)
    arr_loss_delta = np.where(arr_loss_any, arr_loss_delta, 0.0).astype(np.float32)
    arr_loss_idx   = np.where(arr_loss_any, arr_loss_idx, -1).astype(np.int32)
    arr_loss_obs   = np.where(arr_loss_any, arr_loss_obs, np.nan).astype(np.float32)
    arr_loss_ubs   = np.where(arr_loss_any, arr_loss_ubs, np.nan).astype(np.float32)

    dct_out = {
        "Gain_Any":   arr_gain_any,
        "Loss_Any":   arr_loss_any,
        "Gain_Delta": arr_gain_delta,
        "Loss_Delta": arr_loss_delta,
        "Gain_Index": arr_gain_idx,
        "Loss_Index": arr_loss_idx,
        "Gain_Obs":   arr_gain_obs,
        "Gain_Ubs":   arr_gain_ubs,
        "Loss_Obs":   arr_loss_obs,
        "Loss_Ubs":   arr_loss_ubs,
    }

    return dct_out


def decode_pos_strand(num_index, num_strands=2):
    """
    Decode flattened index (0..P*2-1) into (pos, strand).
    strand: 0=Forward, 1=Reverse
    """
    if num_index < 0:
        return None, None
    pos = num_index // num_strands
    s   = num_index %  num_strands
    txt_strand = "Forward" if s == 0 else "Reverse"
    return int(pos), txt_strand


def extract_events_from_reduced(
    dct_red,
    arr_index_seq,
    arr_index_motif,
    txt_chunk_prefix
):
    """
    Build event table from reduced (N,M) arrays.

    Returns
    -------
    dat_event : pandas.DataFrame
        One row per (variant, motif) with >=1 gain and/or >=1 loss event.
    """
    arr_gain_any = dct_red["Gain_Any"]
    arr_loss_any = dct_red["Loss_Any"]

    arr_gain_delta = dct_red["Gain_Delta"]
    arr_loss_delta = dct_red["Loss_Delta"]

    arr_gain_idx = dct_red["Gain_Index"]
    arr_loss_idx = dct_red["Loss_Index"]

    arr_gain_obs = dct_red["Gain_Obs"]
    arr_gain_ubs = dct_red["Gain_Ubs"]

    arr_loss_obs = dct_red["Loss_Obs"]
    arr_loss_ubs = dct_red["Loss_Ubs"]

    N, M = arr_gain_any.shape

    lst_records = []

    ### Gain rows
    idx_nm = np.argwhere(arr_gain_any)
    for n, m in idx_nm:
        txt_variant = str(arr_index_seq[n])
        txt_motif   = str(arr_index_motif[m])

        p, txt_strand = decode_pos_strand(int(arr_gain_idx[n, m]))

        lst_records.append({
            "Chunk_Prefix": txt_chunk_prefix,
            "Variant_ID":   txt_variant,
            "Motif_Name":   txt_motif,
            "Position":     p,
            "Strand":       txt_strand,
            "Delta":        float(arr_gain_delta[n, m]),
            "Obs":          float(arr_gain_obs[n, m]),
            "Ubs":          float(arr_gain_ubs[n, m]),
            "Event_Type":   "Gain",
        })

    ### Loss rows
    idx_nm = np.argwhere(arr_loss_any)
    for n, m in idx_nm:
        txt_variant = str(arr_index_seq[n])
        txt_motif   = str(arr_index_motif[m])

        p, txt_strand = decode_pos_strand(int(arr_loss_idx[n, m]))

        lst_records.append({
            "Chunk_Prefix": txt_chunk_prefix,
            "Variant_ID":   txt_variant,
            "Motif_Name":   txt_motif,
            "Position":     p,
            "Strand":       txt_strand,
            "Delta":        float(arr_loss_delta[n, m]),
            "Obs":          float(arr_loss_obs[n, m]),
            "Ubs":          float(arr_loss_ubs[n, m]),
            "Event_Type":   "Loss",
        })

    dat_event = pd.DataFrame(lst_records)
    return dat_event


def summarize_motif_level_from_any(arr_gain_any, arr_loss_any, lst_motif_name):
    """
    Summarize motif-level gain/loss using (N,M) any-masks.
    """
    M = arr_gain_any.shape[1]
    if len(lst_motif_name) != M:
        raise ValueError(f"Motif name length {len(lst_motif_name)} != M={M}")

    dat_motif = pd.DataFrame({
        "Motif_Name": lst_motif_name,
        "Count_Gain": arr_gain_any.sum(axis=0),
        "Count_Loss": arr_loss_any.sum(axis=0),
    })
    return dat_motif


def summarize_variant_level_from_any(arr_gain_any, arr_loss_any, arr_variant_idx):
    """
    Summarize variant-level gain/loss using (N,M) any-masks.
    """
    N = arr_gain_any.shape[0]
    if len(arr_variant_idx) != N:
        raise ValueError(f"Variant index length {len(arr_variant_idx)} != N={N}")

    dat_variant = pd.DataFrame({
        "Variant_ID": arr_variant_idx,
        "Count_Gain": arr_gain_any.sum(axis=1),
        "Count_Loss": arr_loss_any.sum(axis=1),
    })
    return dat_variant


def main(args):
    # ============================
    # Load scan results
    # ----------------------------
    time_start = time.time()

    ### load data
    print(f"Loading scan results: {args.txt_fpath_scan}")
    obj = np.load(args.txt_fpath_scan, allow_pickle=True)

    arr_motif_scan_ref = obj["Ref"]
    arr_motif_scan_obs = obj["Obs"]
    arr_motif_scan_ubs = obj["Ubs"]
    
    arr_index_seq   = obj["Idx_Sequence"]
    arr_index_motif = obj["Idx_Motif"]

    if arr_motif_scan_obs.shape != arr_motif_scan_ubs.shape:
        raise ValueError("Obs/Ubs shape mismatch.")
        
    if arr_motif_scan_ref.shape != arr_motif_scan_obs.shape:
        raise ValueError("Ref/Obs shape mismatch (expected same scan tensor shape).")
    
    ### calculate dimensions
    N, M, P, S = arr_motif_scan_obs.shape
    print(f"N={N}, M={M}, P={P}, strands={S}")

    if S != 2:
        raise ValueError(f"Expected 2 strands in last dimension, got S={S}")
        
    if len(arr_index_seq) != N:
        raise ValueError(f"Idx_Sequence length {len(arr_index_seq)} != N={N}")
        
    if len(arr_index_motif) != M:
        raise ValueError(f"Idx_Motif length {len(arr_index_motif)} != M={M}")
    
    time_runtime = time.time() - time_start
    print(f"Load and check in {time_runtime:.2f} seconds\n")

    # ============================
    # Load Tbind thresholds
    # ----------------------------
    time_start = time.time()

    ### load motif binding threshold
    print(f"Loading motif Tbind: {args.txt_fpath_tbind}")
    with open(args.txt_fpath_tbind, "rb") as f:
        obj = pickle.load(f)
    
    ### expected structure: {"meta": ..., "tbind": {motif_name: num_Tbind}}
    dct_motif_tbind = obj["tbind"]
    
    ### sanity check: motif order must match scan results
    missing = [m for m in arr_index_motif if m not in dct_motif_tbind]
    if missing:
        raise KeyError(f"{len(missing)} motifs from scan not found in tbind; example: {missing[:3]}")
    
    arr_motif_tbind = np.array(
        [dct_motif_tbind[m] for m in arr_index_motif],
        dtype=np.float32
    )
    
    print(f"Loaded Tbind for {arr_motif_tbind.shape[0]} motifs.")
    
    time_runtime = time.time() - time_start
    print(f"Load and check in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Compute reduced delta + events
    # ----------------------------
    time_start = time.time()

    print("Computing reduced best gain/loss per variant × motif...")
    dct_res = reduce_best_gain_loss(
        arr_motif_scan_obs,
        arr_motif_scan_ubs,
        arr_motif_tbind
    )

    time_runtime = time.time() - time_start
    print(f"Compute reduced delta + best events in {time_runtime:.2f} seconds\n")

    # ============================
    # Extract events + summaries
    # ----------------------------
    time_start = time.time()

    txt_chunk_prefix = os.path.basename(args.txt_fpath_output_prefix)

    print("Extracting gain/loss events (one gain + one loss per variant × motif)...")
    dat_event = extract_events_from_reduced(
        dct_res,
        arr_index_seq,
        arr_index_motif,
        txt_chunk_prefix=txt_chunk_prefix,
    )

    print("Summarizing motif-level gain/loss...")
    dat_summary_motif = summarize_motif_level_from_any(
        dct_res["Gain_Any"],
        dct_res["Loss_Any"],
        arr_index_motif
    )
    
    print("Summarizing variant-level gain/loss...")
    dat_summary_variant = summarize_variant_level_from_any(
        dct_res["Gain_Any"],
        dct_res["Loss_Any"],
        arr_index_seq
    )

    time_runtime = time.time() - time_start
    print(f"Events + summaries in {time_runtime:.2f} seconds\n")

    # ============================
    # Save results
    # ----------------------------
    time_start = time.time()

    txt_fpath_out_results = f"{args.txt_fpath_output_prefix}.delta.results.npz"
    txt_fpath_out_event   = f"{args.txt_fpath_output_prefix}.delta.event.tsv"
    txt_fpath_out_motif   = f"{args.txt_fpath_output_prefix}.delta.summary_motif.tsv"
    txt_fpath_out_variant = f"{args.txt_fpath_output_prefix}.delta.summary_variant.tsv"
    
    print(f"Saving reduced delta NPZ to {txt_fpath_out_results}")
    np.savez(
        txt_fpath_out_results,
        Tbind        = arr_motif_tbind,
        Gain_Any     = dct_res["Gain_Any"],
        Loss_Any     = dct_res["Loss_Any"],
        Gain_Delta   = dct_res["Gain_Delta"],
        Loss_Delta   = dct_res["Loss_Delta"],
        Gain_Index   = dct_res["Gain_Index"],
        Loss_Index   = dct_res["Loss_Index"],
        Gain_Obs     = dct_res["Gain_Obs"],
        Gain_Ubs     = dct_res["Gain_Ubs"],
        Loss_Obs     = dct_res["Loss_Obs"],
        Loss_Ubs     = dct_res["Loss_Ubs"],
        Idx_Sequence = arr_index_seq,
        Idx_Motif    = arr_index_motif,
        P            = np.array([P], dtype=np.int32),
    )

    print(f"Saving events (gain/loss) to {txt_fpath_out_event}")
    dat_event.to_csv(txt_fpath_out_event, sep="\t", index=False)

    print(f"Saving motif summary to {txt_fpath_out_motif}")
    dat_summary_motif.to_csv(txt_fpath_out_motif, sep="\t", index=False)

    print(f"Saving variant summary to {txt_fpath_out_variant}")
    dat_summary_variant.to_csv(txt_fpath_out_variant, sep="\t", index=False)

    time_runtime = time.time() - time_start
    print(f"Save results in {time_runtime:.2f} seconds\n")

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute reduced delta + gain/loss per chunk")
    parser.add_argument("--txt_fpath_scan",  required=True, help="Path to motif scan score (npz file)")
    parser.add_argument("--txt_fpath_tbind", required=True, help="Path to motif Tbind pickle file (output of lods2pmap)")
    parser.add_argument(
        "--txt_fpath_output_prefix",
        required=True,
        help=(
            "Output prefix; script adds "
            "_delta_reduced.npz, _summary_motif.tsv, _summary_variant.tsv, "
            "_event.tsv"
        ),
    )
    args = parser.parse_args()
    main(args)
