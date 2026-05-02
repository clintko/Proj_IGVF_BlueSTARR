#!/usr/bin/env python3
import itertools as it
import argparse
import ushuffle

from fun_fasta import fun_read_fasta

def parse_idx(rid: str):
    chrom, pos, ref, obs, ubs = rid.split(":")
    if not (len(ref) == len(obs) == len(ubs) == 1):
        raise ValueError(f"Non-SNV alleles in id: {rid}")
    return chrom, int(pos), ref, obs, ubs

def shuffle_ref_keep_center(txt_seq_ref, num_var_pos):
    txt_seq_ref = txt_seq_ref.upper()
    txt_seqL = txt_seq_ref[:num_var_pos]
    txt_seqR = txt_seq_ref[num_var_pos+1:]
    txt_var_ref = txt_seq_ref[num_var_pos]

    byt_seqL = txt_seqL.encode("ascii")
    byt_seqR = txt_seqR.encode("ascii")

    shufL = ushuffle.Shuffler(byt_seqL, 2) if len(byt_seqL) >= 2 else None
    shufR = ushuffle.Shuffler(byt_seqR, 2) if len(byt_seqR) >= 2 else None

    txt_seqL_shuffle = (shufL.shuffle() if shufL else byt_seqL).decode("ascii")
    txt_seqR_shuffle = (shufR.shuffle() if shufR else byt_seqR).decode("ascii")

    return txt_seqL_shuffle + txt_var_ref + txt_seqR_shuffle

def main(args):
    ushuffle.set_seed(args.num_seed)

    for rid, seq_ref in it.islice(fun_read_fasta(args.txt_finp), 10):
        print(rid)
        print("Seq(Ref):", seq_ref[:50])
        print("Shuffled:", shuffle_ref_keep_center(seq_ref, args.num_pos0)[:50])
        print()

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--txt_finp", required=True)
    p.add_argument("--num_pos0", type=int, required=True)
    p.add_argument("--num_seed", type=int, required=True)
    main(p.parse_args())
