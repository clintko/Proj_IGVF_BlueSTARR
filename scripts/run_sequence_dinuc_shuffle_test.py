#!/usr/bin/env python3
import argparse
import random
from collections import defaultdict, Counter

from fun_fasta import fun_read_fasta


# ============================================================
# Dinucleotide-preserving shuffle (pure Python; stable)
# ------------------------------------------------------------
def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """
    Dinucleotide-preserving shuffle using an Eulerian trail on a multigraph.
    Pure Python and stable.
    """
    seq = seq.upper()
    n = len(seq)
    if n <= 2:
        return seq

    outgoing = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        outgoing[a].append(b)

    for a in outgoing:
        rng.shuffle(outgoing[a])

    start = seq[0]

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
    if len(path) != n:
        # very rare fallback
        return seq

    return "".join(path)


def shuffle_ref_keep_center(seq_ref: str, pos0: int, rng: random.Random) -> str:
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
# Optional validation helpers
# ------------------------------------------------------------
def dinuc_counts(seq: str) -> Counter:
    seq = seq.upper()
    if len(seq) < 2:
        return Counter()
    return Counter(a + b for a, b in zip(seq[:-1], seq[1:]))


def validate_flanks(seq_ref: str, seq_shuf: str, pos0: int) -> tuple[bool, str]:
    """
    Check dinucleotide counts are identical within left flank and right flank.
    """
    L_ref = seq_ref[:pos0]
    R_ref = seq_ref[pos0 + 1 :]
    L_shu = seq_shuf[:pos0]
    R_shu = seq_shuf[pos0 + 1 :]

    if dinuc_counts(L_ref) != dinuc_counts(L_shu):
        return False, "Left flank dinucleotide counts differ"
    if dinuc_counts(R_ref) != dinuc_counts(R_shu):
        return False, "Right flank dinucleotide counts differ"
    if seq_ref[pos0].upper() != seq_shuf[pos0].upper():
        return False, "Center base changed"
    if len(seq_ref) != len(seq_shuf):
        return False, "Length changed"
    return True, "OK"


# ============================================================
# Main
# ------------------------------------------------------------
def main(args):
    rng = random.Random(args.num_seed)

    n = 0
    for rid, seq_ref in fun_read_fasta(args.txt_finp):
        seq_shuf = shuffle_ref_keep_center(seq_ref, args.num_pos0, rng)

        print(rid)
        print("Seq(Ref): ", seq_ref[: args.num_preview])
        print("Shuffled: ", seq_shuf[: args.num_preview])

        if args.flag_validate:
            ok, msg = validate_flanks(seq_ref, seq_shuf, args.num_pos0)
            print("Validate:", "PASS" if ok else "FAIL", "-", msg)

        print()
        n += 1
        if n >= args.num_n:
            break

    print(f"Printed {n} records.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Test dinucleotide shuffling (pure Python): print first N sequences and optional validation."
    )
    p.add_argument("--txt_finp", required=True, help="Input FASTA (.fa or .fa.gz)")
    p.add_argument("--num_pos0", type=int, required=True, help="0-based variant position in window")
    p.add_argument("--num_seed", type=int, required=True, help="Random seed")
    p.add_argument("--num_n", type=int, default=10, help="Number of records to print")
    p.add_argument("--num_preview", type=int, default=60, help="Preview length to print per sequence")
    p.add_argument("--flag_validate", action="store_true", help="Validate dinuc counts per flank + center fixed")
    main(p.parse_args())
