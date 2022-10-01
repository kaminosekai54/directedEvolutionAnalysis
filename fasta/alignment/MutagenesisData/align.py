#!/usr/bin/env python
from Bio import pairwise2
from sys import argv
from multiprocessing import Pool
from numpy import mean


def read_fasta(fasta_file):
    results = {}
    for l in open(fasta_file):
        if l.startswith(">"):
            name = l.strip()[1:]
            results[name] = ""
        else:
            results[name] += l.strip()
    return results
    

def get_seq(arg):
    "get the alignment to the reference"
    name, cur_seq = arg
    alig = pairwise2.align.globalms(ref_seq, cur_seq.replace("N", "-"), 2, -2, -4, -1, penalize_end_gaps=False)[0]
    nrjs, alig_scr = [], []

    # for alig in aligs:
    al_pos = "".join([cs for az, cs in zip(alig[0], alig[1]) if az != "-"])
    nseq = al_pos + ref_seq
    return name, nseq, alig[-1]


if __name__ == '__main__':
    infile = argv[1]
    h_seq = read_fasta(infile)
    ref_seq = "GTGGTAAAATCTGCCTAAACGGGGAAACTCTCACTGAGACAATCCCGTGCTAAATCAGCAATAGCTGTAAATGCCTAACGACTACACGGTAGACAACTCTAAGAGTTGAAGGTATAGTCTAAACTGCAAGGTGACTTGCAGATATCGG"

    pool = Pool(30)
    seq_list = [(n, seq) for n, seq in h_seq.items()]
    scores = pool.map(get_seq, ((n, s) for n, s in seq_list))
    pool.close()

    for name, seq, ali_scr in scores:
        nb_mut = sum(1 for s, az in zip(seq, ref_seq) if s != az and s != "-")
        print(f">{name} SCORE={ali_scr} NB_MUT={nb_mut}\n{seq}")
