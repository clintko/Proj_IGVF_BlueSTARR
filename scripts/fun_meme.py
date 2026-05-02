import numpy


def read_meme(txt_fpath, num_motifs=None):
    """Read a MEME file and return a dictionary of PWMs.
    This file is modified from memelite

    Parameters
    ----------
    txt_fpath: str
        The filename of the MEME-formatted file to read in

    Returns
    -------
    motifs: dict
        A dictionary of the motifs in the MEME file.
    """

    motifs = {}

    with open(txt_fpath, "r") as file:
        motif, width, i = None, None, 0

        for line in file:
            if motif is None:
                if line[:5] == 'MOTIF':
                    motif = line.replace('MOTIF ', '').strip("\r\n")
                else:
                    continue

            elif width is None:
                if line[:6] == 'letter':
                    width = int(line.split()[5])
                    pwm = numpy.zeros((width, 4))

            elif i < width:
                pwm[i] = list(map(float, line.strip("\r\n").split()))
                i += 1

            else:
                motifs[motif] = pwm #motifs[motif] = pwm.T
                motif, width, i = None, None, 0

                if num_motifs is not None and len(motifs) == num_motifs:
                    break

    return motifs