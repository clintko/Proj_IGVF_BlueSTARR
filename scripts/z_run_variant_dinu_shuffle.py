#!/usr/bin/env python3
import argparse
import ushuffle
import random

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
    txt_seq  = txt_seq_ref.upper().encode("ascii")
    txt_seqL = txt_seq[:num_var_pos]
    txt_seqR = txt_seq[num_var_pos+1:]
    txt_var_ref = txt_seq[num_var_pos:num_var_pos+1]
    
    shuffleL = ushuffle.Shuffler(txt_seqL, 2) if len(txt_seqL)  >= 2 else None
    shuffleR = ushuffle.Shuffler(txt_seqR, 2) if len(txt_seqR) >= 2 else None

    txt_seqL_shuffle = shuffleL.shuffle().decode("ascii")
    txt_seqR_shuffle = shuffleR.shuffle().decode("ascii")

    txt_seq_shuffle_ref = (txt_seqL_shuffle + txt_var_ref + txt_seqR_shuffle)
    txt_seq_shuffle_obs = (txt_seqL_shuffle + txt_var_obs + txt_seqR_shuffle)
    txt_seq_shuffle_ubs = (txt_seqL_shuffle + txt_var_ubs + txt_seqR_shuffle)
    return txt_seq_shuffle_ref, txt_seq_shuffle_obs, txt_seq_shuffle_ubs

def main(args):

    #lst_seq_shuffle_ref = list()
    #lst_seq_shuffle_obs = list()
    #lst_seq_shuffle_ubs = list()
    
    txt_fout_ref = f"{args.txt_fout_prefix}_ref.fa"
    txt_fout_obs = f"{args.txt_fout_prefix}_obs.fa"
    txt_fout_ubs = f"{args.txt_fout_prefix}_ubs.fa"
    
    for txt_seq_idx, txt_seq_ref in read_fasta(args.txt_fpath_ref):
        txt_chrom, num_pos, txt_var_ref, txt_var_obs, txt_var_ubs = parse_idx(txt_seq_idx)
        txt_seq_shuffle_ref, txt_seq_shuffle_obs, txt_seq_shuffle_ubs = shuffle(txt_seq_ref, txt_var_obs, txt_var_ubs, num_var_pos)

        write_fasta(txt_fout_ref, txt_seq_idx, lst_seq_shuffle_ref)
        write_fasta(txt_fout_obs, txt_seq_idx, lst_seq_shuffle_obs)
        write_fasta(txt_fout_ubs, txt_seq_idx, lst_seq_shuffle_ubs)

        #lst_seq_shuffle_ref.append(txt_seq_shuffle_ref)
        #lst_seq_shuffle_obs.append(txt_seq_shuffle_obs)
        #lst_seq_shuffle_ubs.append(txt_seq_shuffle_ubs)
        
        #write_fasta(lst_seq_shuffle_ref, txt_fout_ref)
        #write_fasta(lst_seq_shuffle_obs, txt_fout_obs)
        #write_fasta(lst_seq_shuffle_ubs, txt_fout_ubs)
        

if __name__ == "__main__":
    ### parse arguments
    parser = argparse.ArgumentParser(description="Run motif scanning on variant FASTA files")
    
    parser.add_argument("--txt_finp",        type=str, required=True, help="Path to input fasta file")
    parser.add_argument("--txt_fout_prefix", type=str, required=True, help="Path to output fasta file")
    parser.add_argument("--num_pos0",        type=int, required=True, help="Number of position")
    #parser.add_argument("--num_seed",        type=int, required=True, help="Random seed")

    args = parser.parse_args()
    
    ### run main function
    main(args)