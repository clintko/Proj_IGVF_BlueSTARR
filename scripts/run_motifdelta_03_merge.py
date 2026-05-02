"""
Merge per-chunk motifdelta outputs into pilot-level merged tables.

Expected per-chunk inputs (directory):
    *.delta.event.tsv
    *.delta.summary_motif.tsv
    *.delta.summary_variant.tsv

Outputs (prefix):
    <prefix>.pilot_merge.delta.event.tsv
    <prefix>.pilot_merge.delta.summary_motif.tsv
    <prefix>.pilot_merge.delta.summary_variant.tsv
"""

import pandas as pd
import argparse
import glob
import os
import time

def merge_event_tables(lst_fpath_event):
    """Row-bind all per-chunk event tables."""
    if not lst_fpath_event:
        return pd.DataFrame()

    lst_dat = []
    for fp in lst_fpath_event:
        dat = pd.read_csv(fp, sep="\t")
        if "Delta" in dat.columns:
            dat["Delta"] = pd.to_numeric(dat["Delta"], errors="coerce")
        lst_dat.append(dat)

    return pd.concat(lst_dat, ignore_index=True)


def merge_motif_summaries(lst_fpath_motif):
    """Sum motif-level gain/loss counts across chunks."""
    if not lst_fpath_motif:
        return pd.DataFrame(columns=["Motif_Name", "Count_Gain", "Count_Loss"])

    dat_all = pd.concat(
        [pd.read_csv(fp, sep="\t") for fp in lst_fpath_motif],
        ignore_index=True,
    )

    return (
        dat_all
        .groupby("Motif_Name", as_index=False)[["Count_Gain", "Count_Loss"]]
        .sum()
    )


def merge_variant_summaries(lst_fpath_variant):
    """Sum variant-level gain/loss counts across chunks."""
    if not lst_fpath_variant:
        return pd.DataFrame(columns=["Variant_ID", "Count_Gain", "Count_Loss"])

    dat_all = pd.concat(
        [pd.read_csv(fp, sep="\t") for fp in lst_fpath_variant],
        ignore_index=True,
    )

    return (
        dat_all
        .groupby("Variant_ID", as_index=False)[["Count_Gain", "Count_Loss"]]
        .sum()
    )

#import re
#def list_chunk_files(txt_dir, suffix):
#    fps = sorted(glob.glob(os.path.join(txt_dir, f"*{suffix}")))
#    return [fp for fp in fps if re.search(r"\.chunk\d+\.", os.path.basename(fp))]
def list_chunk_files(txt_dir, suffix):
    fps = sorted(glob.glob(os.path.join(txt_dir, f"*{suffix}")))
    return [fp for fp in fps if "chunk" in os.path.basename(fp)]


def main(args):
    t0 = time.time()
    txt_dir = args.txt_fpath_dir

    #lst_event   = sorted(glob.glob(os.path.join(txt_dir, "*.delta.event.tsv")))
    #lst_motif   = sorted(glob.glob(os.path.join(txt_dir, "*.delta.summary_motif.tsv")))
    #lst_variant = sorted(glob.glob(os.path.join(txt_dir, "*.delta.summary_variant.tsv")))
    lst_event   = list_chunk_files(txt_dir, ".delta.event.tsv")
    lst_motif   = list_chunk_files(txt_dir, ".delta.summary_motif.tsv")
    lst_variant = list_chunk_files(txt_dir, ".delta.summary_variant.tsv")
    
    print(f"Found {len(lst_event)} event files.")
    print(f"Found {len(lst_motif)} motif summary files.")
    print(f"Found {len(lst_variant)} variant summary files.\n")

    print("Merging events...")
    dat_event = merge_event_tables(lst_event)

    print("Merging motif summaries...")
    dat_motif = merge_motif_summaries(lst_motif)

    print("Merging variant summaries...")
    dat_variant = merge_variant_summaries(lst_variant)

    prefix = args.txt_fpath_output_prefix

    fp_out_event   = f"{prefix}.delta.event.tsv"
    fp_out_motif   = f"{prefix}.delta.summary_motif.tsv"
    fp_out_variant = f"{prefix}.delta.summary_variant.tsv"

    print("\nSaving outputs:")
    print(f"  {fp_out_event}")
    dat_event.to_csv(fp_out_event, sep="\t", index=False)

    print(f"  {fp_out_motif}")
    dat_motif.to_csv(fp_out_motif, sep="\t", index=False)

    print(f"  {fp_out_variant}")
    dat_variant.to_csv(fp_out_variant, sep="\t", index=False)

    print(f"\nDone in {time.time() - t0:.2f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge motifdelta pilot chunk outputs (delta naming)."
    )
    parser.add_argument(
        "--txt_fpath_dir",
        required=True,
        help="Directory containing per-chunk *.delta.*.tsv files.",
    )
    parser.add_argument(
        "--txt_fpath_output_prefix",
        required=True,
        help="Output prefix (path + base name, without suffix).",
    )
    args = parser.parse_args()
    main(args)
