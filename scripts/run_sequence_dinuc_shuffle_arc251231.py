#!/usr/bin/env python3
import itertools as it
import argparse
import ushuffle

from fun_fasta import fun_read_fasta

def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)

def write_fasta(txt_fpath, txt_idx, txt_seq, wrap=80):
    with open(txt_fpath, "a") as out:
        out.write(f">{txt_idx}\n")
        for i in range(0, len(txt_seq), wrap):
            out.write(txt_seq[i:i+wrap] + "\n")

def parse_idx(rid: str):
    # rid example: chr4:74487586:A:G:C
    chrom, pos, ref, obs, ubs = rid.split(":")
    if not (len(ref) == len(obs) == len(ubs) == 1):
        raise ValueError(f"Non-SNV alleles in id: {rid}")
    return chrom, int(pos), ref, obs, ubs

def shuffle(txt_seq_ref, txt_var_obs, txt_var_ubs, num_var_pos):
    txt_seq_ref = txt_seq_ref.upper()

    ### keep everything as str
    txt_seqL = txt_seq_ref[:num_var_pos]
    txt_seqR = txt_seq_ref[num_var_pos+1:]
    txt_var_ref = txt_seq_ref[num_var_pos]
    
    ### bytes for ushuffle
    byt_seqL = txt_seqL.encode("ascii")
    byt_seqR = txt_seqR.encode("ascii")

    shuffleL = ushuffle.Shuffler(byt_seqL, 2) if len(byt_seqL) >= 2 else None
    shuffleR = ushuffle.Shuffler(byt_seqR, 2) if len(byt_seqR) >= 2 else None

    byt_seqL_shuffle = shuffleL.shuffle() if shuffleL else byt_seqL
    byt_seqR_shuffle = shuffleR.shuffle() if shuffleR else byt_seqR

    txt_seqL_shuffle = byt_seqL_shuffle.decode("ascii")
    txt_seqR_shuffle = byt_seqR_shuffle.decode("ascii")

    txt_seq_shuffle_ref = (txt_seqL_shuffle + txt_var_ref    + txt_seqR_shuffle)
    txt_seq_shuffle_obs = (txt_seqL_shuffle + txt_var_obs    + txt_seqR_shuffle)
    txt_seq_shuffle_ubs = (txt_seqL_shuffle + txt_var_ubs    + txt_seqR_shuffle)

    return txt_seq_shuffle_ref, txt_seq_shuffle_obs, txt_seq_shuffle_ubs

def main(args):

    ushuffle.set_seed(args.num_seed)

    ### test: first 5 sequences
    for txt_seq_idx, txt_seq_ref in it.islice(read_fasta(args.txt_finp), 10):

        txt_chrom, num_pos, txt_var_ref, txt_var_obs, txt_var_ubs = parse_idx(txt_seq_idx)

        txt_seq_shuffle_ref, txt_seq_shuffle_obs, txt_seq_shuffle_ubs = shuffle(
            txt_seq_ref,
            txt_var_obs,
            txt_var_ubs,
            args.num_pos0
        )

        print(txt_seq_idx)
        print("Seq(Ref):", txt_seq_ref[:50])
        print("Shuffled:", txt_seq_shuffle_ref[:50])
        myfasta.fun_write_fasta(txt_seq_shuffle_ref, args.txt_fout)

if __name__ == "__main__":
    ### parse argument
    parser = argparse.ArgumentParser(description="Dinucleotide shuffling for variants")

    parser.add_argument("--txt_finp", type=str, required=True, help="Path to input fasta file")
    parser.add_argument("--txt_fout", type=str, required=True, help="Path to output fasta file (unused in test)")
    parser.add_argument("--num_pos0", type=int, required=True, help="0-based variant position in window")
    parser.add_argument("--num_seed", type=int, required=True, help="Random seed")

    args = parser.parse_args()

    ### run main function
    main(args)