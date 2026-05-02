"""
Merge per-chunk motifdelta outputs into global tables.

Inputs (directory):
    *_event.tsv
    *_summary_motif.tsv
    *_summary_variant.tsv

Outputs (prefix):
    <prefix>_summary_motif_merged.tsv
    <prefix>_summary_variant_merged.tsv
    <prefix>_events_all.tsv
    <prefix>_events_top_gain.tsv
    <prefix>_events_top_loss.tsv
"""

import numpy as np
import pandas as pd
import argparse
import glob
import os
import time

def merge_event_table(lst_fpath, num_top_global):
    """
    Merge all per-chunk event tables.

    Parameters
    ----------
    lst_fpath : list[str]
        Paths to *_event.tsv files.
    num_top_global : int
        Number of top gain / loss events to keep globally.

    Returns
    -------
    dat_all : pandas.DataFrame
        All gain/loss events from all chunks.
    dat_top_gain : pandas.DataFrame
        Top gain events (by Delta descending).
    dat_top_loss : pandas.DataFrame
        Top loss events (by Delta ascending).
    """
    ### import all tables and store as a list
    lst_dat = []
    for txt_fpath in lst_fpath:
        dat = pd.read_csv(txt_fpath, sep="\t")
        lst_dat.append(dat)

    ### check non-empty
    if len(lst_dat) == 0:
        raise ValueError("No event files found.")

    ### stack all per-chunk dataframes on top of each other
    dat_all = pd.concat(lst_dat, ignore_index=True)

    ### make sure types are correct
    dat_all["Delta"] = dat_all["Delta"].astype(float)

    ### split by event type
    dat_gain = dat_all[dat_all["Event_Type"] == "Gain"]
    dat_loss = dat_all[dat_all["Event_Type"] == "Loss"]

    ### global ranking
    dat_top_gain = dat_gain.sort_values("Delta", ascending=False).head(num_top_global)
    dat_top_loss = dat_loss.sort_values("Delta", ascending=True).head(num_top_global)

    return dat_all, dat_top_gain, dat_top_loss


def merge_motif_summary(lst_fpath):
    """
    Merge motif-level summaries across chunks.

    Each input file must have columns:
        Motif_Name, Count_Gain, Count_Loss
    """
    ### import all tables and store as a list
    lst_dat = []
    for txt_fpath in lst_fpath:
        dat = pd.read_csv(txt_fpath, sep="\t")
        lst_dat.append(dat)

    ### check non-empty; if empty, return an empty dataframe
    if len(lst_dat) == 0:
        return pd.DataFrame(columns=["Motif_Name", "Count_Gain", "Count_Loss"])

    ### stack all per-chunk dataframes on top of each other
    dat_concat = pd.concat(lst_dat, ignore_index=True)

    ### group by Motif_Name and sum counts across chunks
    dat_merged = (
        dat_concat
        .groupby("Motif_Name", as_index=False)[["Count_Gain", "Count_Loss"]]
        .sum()
    )

    return dat_merged


def merge_variant_summary(lst_fpath):
    """
    Merge variant-level summaries across chunks.

    Each input file must have columns:
        Variant_ID, Count_Gain, Count_Loss
    """
    ### import all tables and store as a list
    lst_dat = []
    for txt_fpath in lst_fpath:
        dat = pd.read_csv(txt_fpath, sep="\t")
        lst_dat.append(dat)

    ### check non-empty; if empty, return an empty dataframe
    if len(lst_dat) == 0:
        return pd.DataFrame(columns=["Variant_ID", "Count_Gain", "Count_Loss"])

    ### stack all per-chunk dataframes on top of each other
    dat_all = pd.concat(lst_dat, ignore_index=True)

    ### group by Variant_ID and sum counts (defensive; should be one chunk per variant)
    dat_merged = (
        dat_all
        .groupby("Variant_ID", as_index=False)[["Count_Gain", "Count_Loss"]]
        .sum()
    )

    return dat_merged


def main(args):
    # ============================
    # Import all summary tables
    # ----------------------------
    time_start = time.time()
    
    ### get file paths
    txt_dir = args.txt_fpath_dir

    lst_fpath_motif   = sorted(glob.glob(os.path.join(txt_dir, "*_summary_motif.tsv")))
    lst_fpath_variant = sorted(glob.glob(os.path.join(txt_dir, "*_summary_variant.tsv")))
    lst_fpath_event   = sorted(glob.glob(os.path.join(txt_dir, "*_event.tsv")))

    print(f"Found {len(lst_fpath_motif)} motif summary files.")
    print(f"Found {len(lst_fpath_variant)} variant summary files.")
    print(f"Found {len(lst_fpath_event)} event files.")

    ### sanity check file exists
    if len(lst_fpath_motif) == 0:
        print("WARNING: No *_summary_motif.tsv found.")
    if len(lst_fpath_variant) == 0:
        print("WARNING: No *_summary_variant.tsv found.")
    if len(lst_fpath_event) == 0:
        print("WARNING: No *_event.tsv found.")

    time_runtime = time.time() - time_start
    print(f"Import all tables in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Merge tables
    # ----------------------------
    time_start = time.time()
    
    print("Merging motif-level summaries...")
    dat_motif_merged = merge_motif_summary(lst_fpath_motif)

    print("Merging variant-level summaries...")
    dat_variant_merged = merge_variant_summary(lst_fpath_variant)

    print("Merging event tables into global tables...")
    dat_event_all, dat_event_top_gain, dat_event_top_loss = merge_event_table(
        lst_fpath_event,
        num_top_global=args.num_top_events_global,
    )

    time_runtime = time.time() - time_start
    print(f"Merge tables in {time_runtime:.2f} seconds\n")
    
    # ============================
    # Save outputs
    # ----------------------------
    time_start = time.time()
    
    txt_prefix_out        = args.txt_fpath_output_prefix
    txt_fpath_out_motif   = f"{txt_prefix_out}_summary_motif_merged.tsv"
    txt_fpath_out_variant = f"{txt_prefix_out}_summary_variant_merged.tsv"
    txt_fpath_out_events  = f"{txt_prefix_out}_events_all.tsv"
    txt_fpath_out_gain    = f"{txt_prefix_out}_events_top_gain.tsv"
    txt_fpath_out_loss    = f"{txt_prefix_out}_events_top_loss.tsv"

    print(f"Saving merged motif summary to {txt_fpath_out_motif}")
    dat_motif_merged.to_csv(txt_fpath_out_motif, sep="\t", index=False)

    print(f"Saving merged variant summary to {txt_fpath_out_variant}")
    dat_variant_merged.to_csv(txt_fpath_out_variant, sep="\t", index=False)

    print(f"Saving all events to {txt_fpath_out_events}")
    dat_event_all.to_csv(txt_fpath_out_events, sep="\t", index=False)

    print(f"Saving top gain events to {txt_fpath_out_gain}")
    dat_event_top_gain.to_csv(txt_fpath_out_gain, sep="\t", index=False)

    print(f"Saving top loss events to {txt_fpath_out_loss}")
    dat_event_top_loss.to_csv(txt_fpath_out_loss, sep="\t", index=False)

    time_runtime = time.time() - time_start
    print(f"Save outputs in {time_runtime:.2f} seconds\n")
    
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge motifdelta chunk-level summaries and events."
    )
    parser.add_argument(
        "--txt_fpath_dir",
        required=True,
        help="Directory containing per-chunk *_summary_*.tsv and *_event.tsv files.",
    )
    parser.add_argument(
        "--txt_fpath_output_prefix",
        required=True,
        help=(
            "Output prefix; script writes "
            "_summary_motif_merged.tsv, _summary_variant_merged.tsv, "
            "_events_all.tsv, _events_top_gain.tsv, _events_top_loss.tsv"
        ),
    )
    parser.add_argument(
        "--num_top_events_global",
        type=int,
        default=1000,
        help=(
            "Number of top gain / top loss events to keep globally."
        ),
    )
    args = parser.parse_args()
    main(args)