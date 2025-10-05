"""
Module for reverse compliment sequence and matrix
"""
import numpy as np

def reverse_complement_sequence(txt_seq):
    """
    Return the reverse complement of a DNA sequence string.

    Parameters
    ----------
    txt_seq : str
        DNA sequence (e.g. "ACGTN")

    Returns
    -------
    str
        Reverse complement of input sequence (e.g. "NACGT")
    """
    dct_complement = str.maketrans(
        "ACGTacgt",  ### from these alphabets
        "TGCAtgca"   ### to   those alphabets
    )
    return txt_seq.translate(dct_complement)[::-1] ### flip the order as well
    
def reverse_complement_matrix(arr_seq_NxB, txt_alphabet = "ACGT"):
    """
    Reverse-complement a position×base matrix of shape (N,B) (Default: N,4).
    The input matrix reversed along positions and with columns 
    swapped A <-> T and C <-> G.

    Parameters
    ----------
    arr_seq_NxB : 
        one-hot enocded matrix of a sequence with shape (N,B).
        Default: shape (N, 4)

    Returns
    -------
    np.ndarray
        one-hot enocded matrix of the Reverse-complement sequence with shape (N, B)
        Default: shape (N, 4)
    """
    ### sanity check: alphabet must contain A, C, G, T
    required = set("ACGT")
    if not required.issubset(set(txt_alphabet)):
        raise ValueError(f"Alphabet must contain A,C,G,T. Got: {txt_alphabet}")
        
    ### map base -> column index
    dct_base_index = {base: idx for idx, base in enumerate(txt_alphabet)}
    
    ### define complement mapping
    dct_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    
    ### build index list for complement columns, preserving input order
    lst_base_index_rc = [dct_base_index[dct_complement[base]] for base in txt_alphabet]

    ### flip positions + swap bases
    arr_seq_rc = arr_seq_NxB[::-1, lst_base_index_rc]
    return arr_seq_rc