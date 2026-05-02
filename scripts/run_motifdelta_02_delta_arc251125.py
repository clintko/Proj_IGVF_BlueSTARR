"""
Compute motif delta arrays and motif gain/loss based on Tbind.

Input:
    <chunk>_scan.npz
    motif_model.pkl  (contains Tbind thresholds)

Output:
    <chunk>_delta.npz
    <chunk>_delta_motif_summary.tsv
    <chunk>_delta_variant_summary.tsv
    <chunk>_top_event.tsv
"""

import numpy as np
import pandas as pd
import argparse
import pickle
import time
import os


def compute_delta(arr_obs, arr_ubs):
    """
    Compute delta scores per motif / position / strand.
    Delta = Ubs - Obs
    
    Parameters
    ----------
    arr_obs : np.ndarray
        Observed-allele scores, shape (N, M, P, 2).
    arr_ubs : np.ndarray
        Unobserved-allele scores, shape (N, M, P, 2).

    Returns
    -------
    arr_delta : np.ndarray
        Delta scores, shape (N, M, P, 2).
    """
    ### check dimension
    if arr_obs.shape != arr_ubs.shape:
        raise ValueError(
            f"Shape mismatch between Obs {arr_obs.shape} and Ubs {arr_ubs.shape}"
        )
        
    ### calculate, return delta
    return arr_ubs - arr_obs


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


def summarize_motif_level(arr_gain, arr_loss, lst_motif_name):
    """
    Summarize motif-level gain/loss events.

    Parameters
    ----------
    arr_gain : np.ndarray (bool)
        Gain mask, shape (N, M, P, 2).
    arr_loss : np.ndarray (bool)
        Loss mask, shape (N, M, P, 2).
    lst_motif_name : list[str]
        Motif names, length M.

    Returns
    -------
    dat_motif : pandas.DataFrame
        Columns:
            Motif_Name
            Count_Gain
            Count_Loss
    """
    if arr_gain.shape != arr_loss.shape:
        raise ValueError(
            f"Gain and loss arrays must have same shape, got {arr_gain.shape} vs {arr_loss.shape}"
        )

    N, M, P, S = arr_gain.shape
    if len(lst_motif_name) != M:
        raise ValueError(
            f"Length of motif name list ({len(lst_motif_name)}) "
            f"does not match motif dimension M={M}"
        )

    ### sum across variants, positions, strands
    num_count_gain = arr_gain.sum(axis=(0, 2, 3))   # (M,)
    num_count_loss = arr_loss.sum(axis=(0, 2, 3))   # (M,)

    dat_motif = pd.DataFrame({
        "Motif_Name": lst_motif_name,
        "Count_Gain": num_count_gain,
        "Count_Loss": num_count_loss
    })

    return dat_motif


def summarize_variant_level(arr_gain, arr_loss, lst_variant_idx):
    """
    Summarize variant-level gain/loss events (how many motifs gained/lost).

    Parameters
    ----------
    arr_gain : np.ndarray (bool)
        Gain mask, shape (N, M, P, 2).
    arr_loss : np.ndarray (bool)
        Loss mask, shape (N, M, P, 2).
    lst_variant_idx : array-like
        Variant IDs, length N.

    Returns
    -------
    dat_variant : pandas.DataFrame
        Columns:
            Variant_ID
            Total_Count_Gain
            Total_Count_Loss
    """
    if arr_gain.shape != arr_loss.shape:
        raise ValueError(
            f"Gain and loss arrays must have same shape, got {arr_gain.shape} vs {arr_loss.shape}"
        )

    N = arr_gain.shape[0]
    if len(lst_variant_idx) != N:
        raise ValueError(
            f"Length of variant index ({len(lst_variant_idx)}) "
            f"does not match N={N}"
        )

    ### collapse motif × position × strand → single count per variant
    arr_gain_per_variant = arr_gain.sum(axis=(1, 2, 3))   # (N,)
    arr_loss_per_variant = arr_loss.sum(axis=(1, 2, 3))   # (N,)

    dat_variant = pd.DataFrame({
        "Variant_ID": lst_variant_idx,
        "Count_Gain": arr_gain_per_variant,
        "Count_Loss": arr_loss_per_variant,
    })

    return dat_variant


def extract_gain_loss_events(
    arr_delta,
    arr_gain,
    arr_loss,
    arr_index_seq,
    arr_index_motif,
    txt_chunk_prefix,
    num_top_gain=1000,
    num_top_loss=1000,
):
    """
    Extract per-chunk top gain/loss events.

    Parameters
    ----------
    arr_delta : np.ndarray
        Delta scores, shape (N, M, P, 2).
    arr_gain : np.ndarray (bool)
        Gain mask, shape (N, M, P, 2).
    arr_loss : np.ndarray (bool)
        Loss mask, shape (N, M, P, 2).
    arr_index_seq : array-like
        Variant IDs, length N.
    arr_index_motif : array-like
        Motif names, length M.
    txt_chunk_prefix : str
        Chunk prefix (base name) for this NPZ file.
    num_top_gain : int
        Number of top gain events to keep (per chunk).
    num_top_loss : int
        Number of top loss events to keep (per chunk).

    Returns
    -------
    dat_events : pandas.DataFrame
        Columns:
            Chunk_Prefix
            Variant_ID
            Motif_Name
            Position
            Strand
            Delta
            Event_Type   ("Gain" or "Loss")
    """
    # ================================
    # Get dimensions
    # --------------------------------
    
    ### sanity check dimension of delta scores and gain/loss masks
    if arr_delta.shape != arr_gain.shape or arr_delta.shape != arr_loss.shape:
        raise ValueError(
            f"Delta / Gain / Loss shapes must match, got "
            f"Delta={arr_delta.shape}, Gain={arr_gain.shape}, Loss={arr_loss.shape}"
        )

    ### get dimension parameters
    N, M, P, S = arr_delta.shape

    ### double check dimension parameters
    if len(arr_index_seq) != N:
        raise ValueError(f"Idx_Sequence length {len(arr_index_seq)} != N={N}")
    if len(arr_index_motif) != M:
        raise ValueError(f"Idx_Motif length {len(arr_index_motif)} != M={M}")
    if S != 2:
        raise ValueError(f"Expected S=2 strands, got S={S}")

    ### init output results
    lst_record = []

    # ================================
    # Get gain events
    # --------------------------------
    
    ### get the index and delta of gain events
    ### note that arr_gain is boolean array
    arr_idx_gain   = np.argwhere(arr_gain)   # (#{gain_events}, 4), each row = (n, m, p, s)
    arr_delta_gain = arr_delta[arr_gain]     # (#{gain_events},)

    ### for each event, get variant/motif/strand/delta values
    for (n, m, p, s), delta in zip(arr_idx_gain, arr_delta_gain):
        txt_variant = str(arr_index_seq[n])
        txt_motif   = str(arr_index_motif[m])
        txt_strand  = "Forward" if int(s) == 0 else "Reverse"

        lst_record.append({
            "Chunk_Prefix": txt_chunk_prefix,
            "Variant_ID":   txt_variant,
            "Motif_Name":   txt_motif,
            "Position":     int(p),
            "Strand":       txt_strand,
            "Delta":        float(delta),
            "Event_Type":   "Gain",
        })

    # ================================
    # Get loss events
    # --------------------------------

    ### get the index and delta of loss events
    ### note that arr_loss is boolean array
    arr_idx_loss   = np.argwhere(arr_loss)   # (#loss_events, 4), each row = (n, m, p, s)
    arr_delta_loss = arr_delta[arr_loss]     # (#loss_events,)
    
    ### for each event, get variant/motif/strand/delta values
    for (n, m, p, s), dval in zip(arr_idx_loss, arr_delta_loss):
        txt_variant = str(arr_index_seq[n])
        txt_motif   = str(arr_index_motif[m])
        txt_strand  = "Forward" if int(s) == 0 else "Reverse"

        lst_record.append({
            "Chunk_Prefix": txt_chunk_prefix,
            "Variant_ID":   txt_variant,
            "Motif_Name":   txt_motif,
            "Position":     int(p),
            "Strand":       txt_strand,
            "Delta":        float(dval),
            "Event_Type":   "Loss",
        })

    if len(lst_record) == 0:
        # no events in this chunk
        return pd.DataFrame(columns=[
            "Chunk_Prefix",
            "Variant_ID",
            "Motif_Name",
            "Position",
            "Strand",
            "Delta",
            "Event_Type",
        ])

    # ================================
    # Rank & keep top events per type
    # --------------------------------

    ### create dataframe
    dat_events = pd.DataFrame(lst_record)

    ### for gains: largest positive Delta
    dat_gain = dat_events[dat_events["Event_Type"] == "Gain"]
    ### for losses: most negative Delta (sort ascending)
    dat_loss = dat_events[dat_events["Event_Type"] == "Loss"]

    ### get the top/sorted delta for gain
    if len(dat_gain) > 0:
        ### sort by delta value
        dat_gain = dat_gain.sort_values("Delta", ascending=False)
        ### if #top events specified, get the top events
        if num_top_gain is not None and num_top_gain > 0:
            dat_gain = dat_gain.head(num_top_gain)

    ### get the top/sorted delta for loss
    if len(dat_loss) > 0:
        ### sort by negative delta value
        dat_loss = dat_loss.sort_values("Delta", ascending=True)
        ### if #top events specified, get the top events
        if num_top_loss is not None and num_top_loss > 0:
            dat_loss = dat_loss.head(num_top_loss)

    ### concatenate the top/sorted gain and loss events
    dat_events_top = pd.concat([dat_gain, dat_loss], ignore_index=True)

    return dat_events_top

    
def main(args):
    # ============================
    # Load scan results
    # ----------------------------

    ### load data
    print(f"Loading scan results: {args.txt_fpath_scan}")
    obj = np.load(args.txt_fpath_scan)

    arr_motif_scan_ref = obj["Ref"]
    arr_motif_scan_obs = obj["Obs"]
    arr_motif_scan_ubs = obj["Ubs"]
    
    arr_index_seq   = obj["Idx_Sequence"]
    arr_index_motif = obj["Idx_Motif"]

    ### calculate dimensions
    N, M, P, S = arr_motif_scan_obs.shape
    print(f"N={N}, M={M}, P={P}, strands={S}")

    if S != 2:
        raise ValueError(f"Expected 2 strands in last dimension, got S={S}")

    print()
    # ============================
    # Load Tbind thresholds
    # ----------------------------
    
    ### load data
    print(f"Loading motif model (Tbind): {args.txt_fpath_model}")
    with open(args.txt_fpath_model, "rb") as f:
        dct_motif_model = pickle.load(f)

    ### use the same motif order as in scan results
    lst_motif_name = [str(txt) for txt in arr_index_motif]

    ### sanity check that all motifs exist in the model
    for txt in lst_motif_name:
        if txt not in dct_motif_model:
            raise KeyError(f"Motif '{txt}' from scan results not found in motif model.")

    ### get Tbind threshold
    ### ensure Tbind vector in the same order as arr_index_motif
    arr_motif_Tbind = np.array(
        [dct_motif_model[txt]["num_Tbind"] for txt in arr_index_motif],
        dtype=float
    )
    print(f"Loaded Tbind for {len(arr_motif_Tbind)} motifs.")
    print()
    
    # ============================
    # Compute delta, gain/loss masks
    # ----------------------------
    
    print("Computing Delta = Ubs − Obs...")
    arr_motif_delta = compute_delta(arr_motif_scan_obs, arr_motif_scan_ubs)

    print("Computing motif gain / loss (per variant × motif × position × strand)...")
    arr_motif_gain, arr_motif_loss = compute_gain_loss(
        arr_motif_scan_obs,
        arr_motif_scan_ubs,
        arr_motif_Tbind
    )

    # ============================
    # Extract top gain/loss events
    # ----------------------------
    print("Extracting per-chunk top gain/loss events...")
    txt_chunk_prefix = os.path.basename(args.txt_fpath_output_prefix)

    dat_top_event = extract_gain_loss_events(
        arr_motif_delta,
        arr_motif_gain,
        arr_motif_loss,
        arr_index_seq,
        arr_index_motif,
        txt_chunk_prefix=txt_chunk_prefix,
        num_top_gain=args.num_top_events,
        num_top_loss=args.num_top_events,
    )
    
    # ============================
    # Summaries results at variant-/motif-levels
    # ----------------------------
    print("Summarizing motif-level gain/loss...")
    dat_summary_motif = summarize_motif_level(
        arr_motif_gain,
        arr_motif_loss,
        lst_motif_name
    )

    print("Summarizing variant-level gain/loss...")
    dat_summary_variant = summarize_variant_level(
        arr_motif_gain,
        arr_motif_loss,
        arr_index_seq
    )
    print()
    
    # ============================
    # Save results
    # ----------------------------
    
    ### set output file paths
    txt_fpath_out_delta   = f"{args.txt_fpath_output_prefix}_delta.npz"
    txt_fpath_out_motif   = f"{args.txt_fpath_output_prefix}_summary_motif.tsv"
    txt_fpath_out_variant = f"{args.txt_fpath_output_prefix}_summary_variant.tsv"
    txt_fpath_out_event   = f"{args.txt_fpath_output_prefix}_top_event.tsv"

    ### save results
    print(f"Saving delta NPZ to {txt_fpath_out_delta}\n")
    np.savez_compressed(
        txt_fpath_out_delta,
        Delta        = arr_motif_delta,
        Tbind        = arr_motif_Tbind,
        Gain         = arr_motif_gain,
        Loss         = arr_motif_loss,
        Idx_Sequence = arr_index_seq,
        Idx_Motif    = arr_index_motif,
    )

    print(f"Saving top gain/loss events to {txt_fpath_out_event}\n")
    dat_top_event.to_csv(txt_fpath_out_event, sep="\t", index=False)

    print(f"Saving motif summary to {txt_fpath_out_motif}\n")
    dat_summary_motif.to_csv(txt_fpath_out_motif, sep="\t", index=False)

    print(f"Saving variant summary to {txt_fpath_out_variant}\n")
    dat_summary_variant.to_csv(txt_fpath_out_variant, sep="\t", index=False)
    
    print("Done.")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute delta + gain/loss per chunk")
    parser.add_argument("--txt_fpath_scan",  required=True, help="Path to motif scan score (npz file)")
    parser.add_argument("--txt_fpath_model", required=True, help="Path to motif model (pickle file; contains Tbind)")
    parser.add_argument(
        "--txt_fpath_output_prefix", 
        required=True,
        help=(
            "Output prefix; script adds "
            "_delta.npz, _summary_motif.tsv, _summary_variant.tsv, _top_event.tsv"
        ),
    )
    parser.add_argument(
        "--num_top_events", 
        type=int,
        default=1000,
        help="Number of top gain/loss events to keep per chunk (per type).",
    )
    args = parser.parse_args()
    main(args)

