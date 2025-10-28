### set environment
import numpy as np
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

from motifdelta.scan   import scan_motif_both_strands
from motifdelta.encode import one_hot_encode
from motifdelta.rc     import reverse_complement_matrix

# ===================================
# Base implementation (pure Python)
# -----------------------------------

def scan_motif_base(txt_seq, arr_motif):
    """
    Baseline motif scanning using explicit for-loop.
    """
    X = one_hot_encode(txt_seq)
    L, W = X.shape[0], arr_motif.shape[0]
    scores = np.empty(L - W + 1, dtype=float)
    
    for idx in range(L - W + 1):
        window = X[idx:idx + W, :]    # <-- fix: use idx, not i
        scores[idx] = np.sum(window * arr_motif)
    return scores


def scan_motif_both_strands_base(txt_seq, arr_motif_WxB, txt_alphabet="ACGT"):
    """
    Scan motif on both forward and reverse-complement strands (base version).
    """
    arr_score_fwd = scan_motif_base(txt_seq, arr_motif_WxB)
    arr_motif_rc  = reverse_complement_matrix(arr_motif_WxB, txt_alphabet=txt_alphabet)
    arr_score_rev = scan_motif_base(txt_seq, arr_motif_rc)
    return arr_score_fwd, arr_score_rev


def scan_one_sequence_base(seq_record, dct_arr_motif_Wx4):
    """
    Baseline motif scanning for one sequence using the pure loop version.
    """
    txt_seq_id  = seq_record.id
    txt_seq_str = str(seq_record.seq)
    dct_results = {}

    for txt_motif_name, arr_motif_Wx4 in dct_arr_motif_Wx4.items():
        arr_scores_fwd, arr_scores_rev = scan_motif_both_strands_base(txt_seq_str, arr_motif_Wx4)
        dct_results[txt_motif_name] = {
            "forward": arr_scores_fwd,
            "reverse": arr_scores_rev
        }
    return txt_seq_id, dct_results

# ===================================
# NumPy-optimized implementation
# -----------------------------------

def scan_one_sequence(seq_record, dct_arr_motif_Wx4):
    """
    Scan all motifs on a single sequence.
    Returns (seq_id, {motif_name: {"forward": fwd_scores, "reverse": rev_scores}}).

    Parameters
    ----------
    seq_record : SeqRecord or object with `.id` and `.seq`
        Sequence record, each with `.id` and `.seq` attributes.
    dct_arr_motif_Wx4 : dict
        Dictionary of motif_name -> motif matrix (W,4)
        
    Returns
    -------
    tuple
        (seq_id, {motif_name: {"forward": arr, "reverse": arr}})
    """
    ### init: parse sequence id and string
    txt_seq_idx = seq_record.id
    txt_seq_str = str(seq_record.seq)
    
    ### init: collect results
    dct_results = {}

    ### loop through each motif and scan
    for txt_motif_name, arr_motif_Wx4 in dct_arr_motif_Wx4.items():
        arr_scores_fwd, arr_scores_rev = scan_motif_both_strands(txt_seq_str, arr_motif_Wx4)
        dct_results[txt_motif_name] = {
            "forward": arr_scores_fwd,
            "reverse": arr_scores_rev
        }

    return txt_seq_idx, dct_results


# ===================================
# Benchmark three batch scanning cases
#     Base:     pure Python base version
#     Serial:   NumPy-optimized, single-core version
#     Parallel: NumPy-optimized + multi-core version
# -----------------------------------

def scan_sequence_batch_base(lst_seq_record, dct_arr_motif_Wx4):
    """
    Scan a batch of sequences (pure Python base version).
    
    Parameters
    ----------
    lst_seq_record : list[SeqRecord] or list[obj]
        List of sequence records, each with `.id` and `.seq`.
    dct_arr_motif_Wx4 : dict
        Dictionary of motif_name -> motif matrix (W,4)
        
    Returns
    -------
    dict
        {seq_id: {motif_name: {"forward": arr, "reverse": arr}}}
    """
    ### init: final results
    dct_results = {}

    ### loop through each sequence and scan for all motifs in dictionary
    for seq_record in lst_seq_record:
        idx, res = scan_one_sequence_base(seq_record, dct_arr_motif_Wx4)
        dct_results[idx] = res
        
    return dct_results

    
def scan_sequence_batch_serial(lst_seq_record, dct_arr_motif_Wx4):
    """
    Scan a batch of sequences (serial version).
    
    Parameters
    ----------
    lst_seq_record : list[SeqRecord] or list[obj]
        List of sequence records, each with `.id` and `.seq`.
    dct_arr_motif_Wx4 : dict
        Dictionary of motif_name -> motif matrix (W,4)
        
    Returns
    -------
    dict
        {seq_id: {motif_name: {"forward": arr, "reverse": arr}}}
    """
    ### init: final results
    dct_results = {}

    ### loop through each sequence and scan for all motifs in dictionary
    for seq_record in lst_seq_record:
        idx, res = scan_one_sequence(seq_record, dct_arr_motif_Wx4)
        dct_results[idx] = res
        
    return dct_results


def scan_sequence_batch_parallel(lst_seq_record, dct_arr_motif_Wx4, num_workers=None):
    """
    Parallel version of sequence scanning using multiple CPU cores.

    Parameters
    ----------
    lst_seq_record : list[SeqRecord]
        List of sequence records, each with `.id` and `.seq`.
    dct_arr_motif_Wx4 : dict
        Dictionary of motif_name -> motif matrix (W,4)
    num_workers : int, optional
        Number of parallel processes. Defaults to available CPU cores.

    Returns
    -------
    dict
        {seq_id: {motif_name: {"forward": arr, "reverse": arr}}}
    """
    ### init: set number of core for parallelization
    ### default: #{available CPU} - 1
    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 1)

    ### init: final results
    dct_results = {}

    ### loop through each sequence and scan for all motifs in dictionary
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        ### assign workers
        futures = [
            executor.submit(scan_one_sequence, seq_record, dct_arr_motif_Wx4)
            for seq_record in lst_seq_record
        ]
        
        ### execute workers
        for f in as_completed(futures):
            idx, res = f.result()
            dct_results[idx] = res

    return dct_results


