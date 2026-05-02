#!/usr/bin/env python3
import argparse
import random
from collections import defaultdict

from fun_fasta import fun_read_fasta, fun_open_text, fun_wrap_fasta


# ============================================================
# Dinucleotide-preserving shuffle (pure Python; stable)
# ------------------------------------------------------------
def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """
    Dinucleotide-preserving shuffle using an Eulerian trail on a multigraph.
    Pure Python and stable (no segfaults).
    """
    seq = seq.upper()
    n = len(seq)
    if n <= 2:
        return seq

    # Build adjacency lists: edges are dinucleotides X->Y
    outgoing = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        outgoing[a].append(b)

    # Randomize edge order per node
    for a in outgoing:
        rng.shuffle(outgoing[a])

    # Start at original first base (keeps mono-nuc composition and dinucs)
    start = seq[0]

    # Hierholzer algorithm for Eulerian path
    stack = [start]
    path = []
    while stack:
        v = stack[-1]
        if outgoing[v]:
            u = outgoing[v].pop()
            stack.append(u)
        else:
            path.append(stack.pop())

    path.reverse()

    # If something went weird, fall back to original (should be rare)
    if len(path) != n:
        return seq

    return "".join(path)


def shuffle_ref_keep_center(seq_ref: str, pos0: int, rng: random.Random) -> str:
    """
    Shuffle left and right flanks with dinucleotide-preserving shuffle separately.
    Keep the center base fixed.
    """
    seq_ref = seq_ref.upper()
    if pos0 < 0 or pos0 >= len(seq_ref):
        raise ValueError(f"num_pos0 out of range: pos0={pos0}, len={len(seq_ref)}")

    left = seq_ref[:pos0]
    center = seq_ref[pos0]
    right = seq_ref[pos0 + 1 :]

    left_s = dinuc_shuffle(left, rng) if len(left) >= 2 else left
    right_s = dinuc_shuffle(right, rng) if len(right) >= 2 else right

    return left_s + center + right_s


# ============================================================
# Main
# ------------------------------------------------------------
def main(args):
    rng = random.Random(args.num_seed)

    n = 0
    with fun_open_text(args.txt_fout, mode="wt") as fout:
        for rid, seq_ref in fun_read_fasta(args.txt_finp):
            seq_shuf = shuffle_ref_keep_center(seq_ref, args.num_pos0, rng)

            fout.write(f">{rid}\n")
            fout.write(fun_wrap_fasta(seq_shuf, num_width=args.num_wrap))
            fout.write("\n")

            n += 1
            if args.num_test > 0 and n >= args.num_test:
                break
            if args.num_verbose_every > 0 and (n % args.num_verbose_every) == 0:
                print(f"Processed {n} sequences...")

    print(f"Done. Wrote {n} records to: {args.txt_fout}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Dinucleotide shuffling (REF only). Shuffle flanks separately; keep center base fixed."
    )
    p.add_argument("--txt_finp", type=str, required=True, help="Input FASTA (.fa or .fa.gz)")
    p.add_argument("--txt_fout", type=str, required=True, help="Output shuffled REF FASTA (.fa or .fa.gz)")
    p.add_argument("--num_pos0", type=int, required=True, help="0-based variant position in window")
    p.add_argument("--num_seed", type=int, required=True, help="Random seed")
    p.add_argument("--num_test", type=int, default=0, help="If >0, stop after N records (debug)")
    p.add_argument("--num_verbose_every", type=int, default=10000, help="Progress print frequency (0 disables)")
    p.add_argument("--num_wrap", type=int, default=60, help="FASTA line wrap width")
    main(p.parse_args())
