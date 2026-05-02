import gzip

def fun_open_text(txt_fpath, mode="rt"):
    """
    Open a text file that may be plain or gzipped.
    Mode: "rt" for read-text, "wt" for write-text, "at" for append-text
    """
    if txt_fpath.endswith(".gz"):
        return gzip.open(txt_fpath, mode)
    else:
        return open(txt_fpath, mode)

def fun_wrap_fasta(txt_seq, num_width=60):
    return "\n".join(
        txt_seq[i:i+num_width]
        for i in range(0, len(txt_seq), num_width)
    )

def fun_read_fasta(txt_fpath):
    """
    Yield (header, seq) from a FASTA file.
    - header: without the leading '>'
    - seq: concatenated sequence string (no newlines)
    """
    ### init
    header = None
    seq_chunks = []

    ### loop through each line of the file
    with fun_open_text(txt_fpath, mode="rt") as file:
        for line in file:
            ### trim out newline and spaces
            line = line.strip()

            ### if the line is empty, skip to the next line
            if not line:
                continue

            ### check and get the fasta header
            if line.startswith(">"):

                ### if we already had a previous header:
                ###   that means we just finished reading its sequence
                if header is not None:
                    yield header, "".join(seq_chunks)

                ### removes the leading `>` to store only the identifier
                ### this becomes the new current header
                header = line[1:]  # drop '>'

                ### resets the sequence accumulator
                seq_chunks = []
            else:
                ### if sequence line (not header), append it to the list of sequence chunks
                seq_chunks.append(line)
                
        ### output the final fasta record at the end of file
        if header is not None:
            yield header, "".join(seq_chunks)

def fun_write_fasta(txt_fpath, records, width=60):
    """
    Write FASTA records to file.

    Parameters
    ----------
    txt_fpath : str
        Output FASTA path
    records : iterable of (header, seq)
        FASTA records
    width : int
        Line width for wrapping
    """
    with fun_open_text(txt_fpath, mode="wt") as fout:
        for header, seq in records:
            fout.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fout.write(seq[i:i+width] + "\n")
