"""
Compute motif delta arrays and motif gain/loss based on Tbind.

Input:
    <chunk>_scan.npz
    motif_model.pkl  (contains Tbind thresholds)

Output:
    <chunk>_delta.npz
    <chunk>_summary_motif.tsv
    <chunk>_summary_variant.tsv
    <chunk>_event.tsv  (combined gain + loss events)
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

    For each (variant, motif) pair, we first ask whether there is
    any gain/loss event at any position/strand. Then we count, per motif,
    how many variants show gain / loss (no double counting per motif–variant).

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
            Count_Gain (number of variants with >=1 gain for this motif)
            Count_Loss (number of variants with >=1 loss for this motif)
    """
    ### sanity check: gain and loss mask should have the same shape
    if arr_gain.shape != arr_loss.shape:
        raise ValueError(
            f"Gain and loss arrays must have same shape, got {arr_gain.shape} vs {arr_loss.shape}"
        )

    ### get the dimensions
    N, M, P, S = arr_gain.shape
    if len(lst_motif_name) != M:
        raise ValueError(
            f"Length of motif name list ({len(lst_motif_name)}) "
            f"does not match motif dimension M={M}"
        )

    ### for each motif-variant pair, ask if there is a gain/loss event
    arr_gain_any = arr_gain.any(axis=(2, 3))  # (N, M) bool
    arr_loss_any = arr_loss.any(axis=(2, 3))  # (N, M) bool

    ### sum across variants, positions, strands
    arr_gain_per_motif = arr_gain_any.sum(axis=0)   # (M,)
    arr_loss_per_motif = arr_loss_any.sum(axis=0)   # (M,)

    ### arrange results into a table
    dat_motif = pd.DataFrame({
        "Motif_Name": lst_motif_name,
        "Count_Gain": arr_gain_per_motif,
        "Count_Loss": arr_loss_per_motif
    })

    return dat_motif


def summarize_variant_level(arr_gain, arr_loss, lst_variant_idx):
    """
    Summarize variant-level gain/loss events.

    For each (variant, motif) pair, we first ask whether there is
    any gain/loss event at any position/strand. Then we count, per variant,
    how many motifs show gain / loss (no double counting per motif).

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
            Count_Gain (number of motifs with >=1 gain in this variant)
            Count_Loss (number of motifs with >=1 loss in this variant)
    """
    ### sanity check: gain and loss mask should have the same shape
    if arr_gain.shape != arr_loss.shape:
        raise ValueError(
            f"Gain and loss arrays must have same shape, got {arr_gain.shape} vs {arr_loss.shape}"
        )

    ### get the number of variants
    N = arr_gain.shape[0]
    if len(lst_variant_idx) != N:
        raise ValueError(
            f"Length of variant index ({len(lst_variant_idx)}) "
            f"does not match N={N}"
        )

    ### for each motif-variant pair, ask if there is a gain/loss event
    arr_gain_any = arr_gain.any(axis=(2,3)) # (N,M) bool
    arr_loss_any = arr_loss.any(axis=(2,3)) # (N,M) bool

    ### collapse motif × position × strand
    ### sum across motif and get single count per variant
    arr_gain_per_variant = arr_gain_any.sum(axis=1)   # (N,)
    arr_loss_per_variant = arr_loss_any.sum(axis=1)   # (N,)

    ### arrange results into a table
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
    txt_chunk_prefix
):
    """
    Extract gain/loss events for each (variant, motif) pair, keeping at most
    one gain (max delta) and one loss (min delta) per pair.

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

    Returns
    -------
    dat_event : pandas.DataFrame
        One row per (variant, motif) with >=1 gain and/or >=1 loss event.
        Columns:
            Chunk_Prefix
            Variant_ID
            Motif_Name
            Position
            Strand
            Delta
            Event_Type  ("Gain" or "Loss")
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
    lst_gain_records = []
    lst_loss_records = []

    # ================================
    # Gain events: pick max delta per (variant, motif)
    # --------------------------------
    
    ### get the delta of gain events; -inf where no gain
    arr_delta_gain_all = np.where(arr_gain, arr_delta, -np.inf)  # (N,M,P,S)

    ### get the max delta and its index for each variant-motif pair
    arr_delta_gain_max = arr_delta_gain_all.max(axis=(2,3)) # (N,M)
    arr_index_gain_max = arr_delta_gain_all.reshape(N, M, -1).argmax(axis=2)  # (N,M)

    ### where there is at least one gain event
    arr_gain_any = arr_gain.any(axis=(2,3)) # (N,M) bool

    ### for each event, get variant/motif/strand/delta values
    for n in range(N):
        txt_variant = str(arr_index_seq[n])
        for m in range(M):
            if not arr_gain_any[n, m]:
                continue   # no gain events for this variant-motif

            delta = arr_delta_gain_max[n, m]
            if not np.isfinite(delta):
                continue   # defensive; should not happen if gain_any is True

            idx_max = int(arr_index_gain_max[n, m])   # 0..(P*S-1)
            p = idx_max // S                          # position 0..P-1
            s = idx_max % S                           # strand 0 or 1

            txt_motif  = str(arr_index_motif[m])
            txt_strand = "Forward" if s == 0 else "Reverse"

            lst_gain_records.append({
                "Chunk_Prefix": txt_chunk_prefix,
                "Variant_ID":   txt_variant,
                "Motif_Name":   txt_motif,
                "Position":     int(p),
                "Strand":       txt_strand,
                "Delta":        float(delta),
                "Event_Type":   "Gain",
            })
                                                
    # ================================
    # Loss events: pick min delta per (variant, motif)
    # --------------------------------
    ### get the delta of loss events; -inf where no loss
    arr_delta_loss_all = np.where(arr_loss, arr_delta, np.inf)  # (N,M,P,S)

    ### get the min delta and its index (over P*S) for each variant-motif pair
    arr_delta_loss_min = arr_delta_loss_all.min(axis=(2, 3)) # (N,M)
    arr_index_loss_min = arr_delta_loss_all.reshape(N, M, -1).argmin(axis=2) # (N,M)

    ### where there is at least one loss event
    arr_loss_any = arr_loss.any(axis=(2,3)) # (N,M) bool
    
    ### for each event, get variant/motif/strand/delta values
    for n in range(N):
        txt_variant = str(arr_index_seq[n])
        for m in range(M):
            if not arr_loss_any[n, m]:
                continue   # no loss events for this variant-motif

            delta = arr_delta_loss_min[n, m]
            if not np.isfinite(delta):
                continue   # defensive; should not happen if loss_any is True

            idx_min = int(arr_index_loss_min[n, m])  # 0..(P*S-1)
            p = idx_min // S                         # position 0..P-1
            s = idx_min % S                          # strand 0 or 1

            txt_motif  = str(arr_index_motif[m])
            txt_strand = "Forward" if s == 0 else "Reverse"

            lst_loss_records.append({
                "Chunk_Prefix": txt_chunk_prefix,
                "Variant_ID":   txt_variant,
                "Motif_Name":   txt_motif,
                "Position":     int(p),
                "Strand":       txt_strand,
                "Delta":        float(delta),
                "Event_Type":   "Loss",
            })

    # ================================
    # combine events into a dataframe
    # --------------------------------

    ### create dataframe
    lst_event = lst_gain_records + lst_loss_records
    dat_event = pd.DataFrame(lst_event)
    
    return dat_event

    
def main(args):
    # ============================
    # Load scan results
    # ----------------------------
    time_start = time.time()
    
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

    time_runtime = time.time() - time_start
    print(f"Load and check in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Load Tbind thresholds
    # ----------------------------
    time_start = time.time()
    
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
    
    time_runtime = time.time() - time_start
    print(f"Load and check in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Compute delta, gain/loss masks
    # ----------------------------
    time_start = time.time()
    
    print("Computing Delta = Ubs − Obs...")
    arr_motif_delta = compute_delta(arr_motif_scan_obs, arr_motif_scan_ubs)

    print("Computing motif gain / loss (per variant × motif × position × strand)...")
    arr_motif_gain, arr_motif_loss = compute_gain_loss(
        arr_motif_scan_obs,
        arr_motif_scan_ubs,
        arr_motif_Tbind
    )

    time_runtime = time.time() - time_start
    print(f"Compute delta, gain/loss masks in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Extract per-(variant, motif) max gain / min loss
    # ----------------------------
    time_start = time.time()
    
    print("Extracting gain/loss events (max gain and min loss per variant × motif)...")
    txt_chunk_prefix = os.path.basename(args.txt_fpath_output_prefix)

    dat_event = extract_gain_loss_events(
        arr_motif_delta,
        arr_motif_gain,
        arr_motif_loss,
        arr_index_seq,
        arr_index_motif,
        txt_chunk_prefix=txt_chunk_prefix,
    )

    time_runtime = time.time() - time_start
    print(f"Extract per-(variant, motif) gain/loss in {time_runtime:.2f} seconds\n")

    # ============================
    # Summaries at motif-/variant-level
    # ----------------------------
    time_start = time.time()
    
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
    
    time_runtime = time.time() - time_start
    print(f"Summaries at motif-/variant-level in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Save results
    # ----------------------------
    time_start = time.time()
    
    txt_fpath_out_delta   = f"{args.txt_fpath_output_prefix}_delta.npz"
    txt_fpath_out_motif   = f"{args.txt_fpath_output_prefix}_summary_motif.tsv"
    txt_fpath_out_variant = f"{args.txt_fpath_output_prefix}_summary_variant.tsv"
    txt_fpath_out_event   = f"{args.txt_fpath_output_prefix}_event.tsv"

    print(f"Saving delta NPZ to {txt_fpath_out_delta}")
    np.savez(
        txt_fpath_out_delta,
        Delta        = arr_motif_delta,
        Tbind        = arr_motif_Tbind,
        Gain         = arr_motif_gain,
        Loss         = arr_motif_loss,
        Idx_Sequence = arr_index_seq,
        Idx_Motif    = arr_index_motif,
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
    parser = argparse.ArgumentParser(description="Compute delta + gain/loss per chunk")
    parser.add_argument("--txt_fpath_scan",  required=True, help="Path to motif scan score (npz file)")
    parser.add_argument("--txt_fpath_model", required=True, help="Path to motif model (pickle file; contains Tbind)")
    parser.add_argument(
        "--txt_fpath_output_prefix", 
        required=True,
        help=(
            "Output prefix; script adds "
            "_delta.npz, _summary_motif.tsv, _summary_variant.tsv, "
            "_event.tsv"
        ),
    )
    args = parser.parse_args()
    main(args)


