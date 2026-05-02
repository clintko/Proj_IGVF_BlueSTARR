"""
Merge motifdelta chunk-level outputs (motif/variant summaries + top events).

Input directory should contain per-chunk files:
    <prefix>_summary_motif.tsv
    <prefix>_summary_variant.tsv
    <prefix>_top_event.tsv

Output:
    <out_prefix>_summary_motif_merged.tsv
    <out_prefix>_summary_variant_merged.tsv
    <out_prefix>_events_all.tsv
    <out_prefix>_events_top_gain.tsv
    <out_prefix>_events_top_loss.tsv
"""

import numpy as np
import pandas as pd
import argparse
import glob
import os


def merge_motif_summary(lst_fpath_motif):
    """
    Merge motif-level summaries across chunks.

    Each input file must have columns:
        Motif_Name, Count_Gain, Count_Loss
    """
    ### import all tables and store as a list
    lst_dat = []
    for txt_fpath in lst_fpath_motif:
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


def merge_variant_summary(lst_fpath_variant):
    """
    Merge variant-level summaries across chunks.

    Each input file must have columns:
        Variant_ID, Total_Count_Gain, Total_Count_Loss
    """
    ### import all tables and store as a list
    lst_dat = []
    for txt_fpath in lst_fpath_variant:
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


def merge_events(lst_fpath_event, num_top_global=None):
    """
    Merge per-chunk top_event tables and produce global rankings.

    Each input file must have columns:
        Chunk_Prefix, Variant_ID, Motif_Name,
        Position, Strand, Delta, Event_Type ("Gain" or "Loss")

    Parameters
    ----------
    lst_fpath_event : list[str]
        List of *_top_event.tsv file paths.
    num_top_global : int or None
        If not None, keep only top-N gains / losses globally.
        If None, keep all.

    Returns
    -------
    dat_concat : pandas.DataFrame
        All concatenated events.
    dat_gain : pandas.DataFrame
        Global top gain events.
    dat_loss : pandas.DataFrame
        Global top loss events.
    """
    ### import all tables and store as a list
    lst_dat = []
    for txt_fpath in lst_fpath_event:
        dat = pd.read_csv(txt_fpath, sep="\t")
        lst_dat.append(dat)

    ### check non-empty; if empty, return an empty dataframe
    if len(lst_dat) == 0:
        empty_cols = [
            "Chunk_Prefix",
            "Variant_ID",
            "Motif_Name",
            "Position",
            "Strand",
            "Delta",
            "Event_Type",
        ]
        dat_empty = pd.DataFrame(columns=empty_cols)
        return dat_empty, dat_empty, dat_empty

    ### stack all per-chunk dataframes on top of each other
    dat_concat = pd.concat(lst_dat, ignore_index=True)

    ### ensure column types
    if "Position" in dat_concat.columns:
        dat_concat["Position"] = dat_concat["Position"].astype(int)
    if "Delta" in dat_concat.columns:
        dat_concat["Delta"] = dat_concat["Delta"].astype(float)

    ### split by event type (gain and loss)
    dat_gain = dat_concat[dat_concat["Event_Type"] == "Gain"].copy()
    dat_loss = dat_concat[dat_concat["Event_Type"] == "Loss"].copy()

    ### global ranking for gain events
    if len(dat_gain) > 0:
        dat_gain = dat_gain.sort_values("Delta", ascending=False)
        if num_top_global is not None and num_top_global > 0:
            dat_gain = dat_gain.head(num_top_global)

    ### global ranking for loss events
    if len(dat_loss) > 0:
        dat_loss = dat_loss.sort_values("Delta", ascending=True)
        if num_top_global is not None and num_top_global > 0:
            dat_loss = dat_loss.head(num_top_global)

    return dat_concat, dat_gain, dat_loss


def main(args):
    # ============================
    # Collect input files
    # ----------------------------
    
    ### get file paths
    txt_dir = args.txt_fpath_dir

    lst_fpath_motif   = sorted(glob.glob(os.path.join(txt_dir, "*_summary_motif.tsv")))
    lst_fpath_variant = sorted(glob.glob(os.path.join(txt_dir, "*_summary_variant.tsv")))
    lst_fpath_event   = sorted(glob.glob(os.path.join(txt_dir, "*_top_event.tsv")))

    print(f"Found {len(lst_fpath_motif)} motif summary files.")
    print(f"Found {len(lst_fpath_variant)} variant summary files.")
    print(f"Found {len(lst_fpath_event)} event files.")

    ### sanity check file exists
    if len(lst_fpath_motif) == 0:
        print("WARNING: No *_summary_motif.tsv found.")
    if len(lst_fpath_variant) == 0:
        print("WARNING: No *_summary_variant.tsv found.")
    if len(lst_fpath_event) == 0:
        print("WARNING: No *_top_event.tsv found.")

    # ============================
    # Merge tables
    # ----------------------------
    print()
    print("Merging motif-level summaries...")
    dat_motif_merged = merge_motif_summary(lst_fpath_motif)

    print("Merging variant-level summaries...")
    dat_variant_merged = merge_variant_summary(lst_fpath_variant)

    print("Merging top event tables (per-chunk) into global tables...")
    dat_event_all, dat_event_top_gain, dat_event_top_loss = merge_events(
        lst_fpath_event,
        num_top_global=args.num_top_events_global,
    )

    # ============================
    # Save outputs
    # ----------------------------
    txt_prefix_out = args.txt_fpath_output_prefix

    txt_fpath_out_motif   = f"{txt_prefix_out}_summary_motif_merged.tsv"
    txt_fpath_out_variant = f"{txt_prefix_out}_summary_variant_merged.tsv"
    txt_fpath_out_events  = f"{txt_prefix_out}_events_all.tsv"
    txt_fpath_out_gain    = f"{txt_prefix_out}_events_top_gain.tsv"
    txt_fpath_out_loss    = f"{txt_prefix_out}_events_top_loss.tsv"

    print(f"\nSaving merged motif summary to {txt_fpath_out_motif}")
    dat_motif_merged.to_csv(txt_fpath_out_motif, sep="\t", index=False)

    print(f"\nSaving merged variant summary to {txt_fpath_out_variant}")
    dat_variant_merged.to_csv(txt_fpath_out_variant, sep="\t", index=False)

    print(f"\nSaving all events to {txt_fpath_out_events}")
    dat_event_all.to_csv(txt_fpath_out_events, sep="\t", index=False)

    print(f"\nSaving global top gain events to {txt_fpath_out_gain}")
    dat_event_top_gain.to_csv(txt_fpath_out_gain, sep="\t", index=False)

    print(f"\nSaving global top loss events to {txt_fpath_out_loss}")
    dat_event_top_loss.to_csv(txt_fpath_out_loss, sep="\t", index=False)
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge motifdelta chunk-level summaries and top events."
    )
    parser.add_argument(
        "--txt_fpath_dir",
        required=True,
        help="Directory containing per-chunk *_summary_*.tsv and *_top_event.tsv files.",
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
            "Number of top gain / top loss events to keep globally. "
            "Use a large number if you want almost all events."
        ),
    )
    args = parser.parse_args()
    main(args)
