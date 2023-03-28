#!/usr/bin/env python
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sys import argv
from multiprocessing import Pool
from numpy import mean
import time, os


def read_fasta(fasta_file):
    results = {}
    for record in SeqIO.parse(fasta_file, "fasta-2line"):
        results[record.id] = str(record.seq)

    return results
    

def get_seq(arg):
    "get the alignment to the reference"
    name, cur_seq = arg
    alig = pairwise2.align.globalms(ref_seq, cur_seq.replace("N", "-"), 2, -2, -4, -1, penalize_end_gaps=False)[0]
    nrjs, alig_scr = [], []

    # for alig in aligs:
    al_pos = "".join([cs for az, cs in zip(alig[0], alig[1]) ])
    nseq = al_pos 
    return name, nseq, alig[-1]


if __name__ == '__main__':
    infile = argv[1]
    h_seq = read_fasta(infile)
    ref_seq = "GTGGTAAAATCTGCCTAAACGGGGAAACTCTCACTGAGACAATCCCGTGCTAAATCAGCAATAGCTGTAAATGCCTAACGACTACACGGTAGACAACTCTAAGAGTTGAAGGTATAGTCTAAACTGCAAGGTGACTTGCAGATATCGG"
    print("starting to align")
    t1 = time.time()

    pool = Pool(100)
    seq_list = [(n, seq) for n, seq in h_seq.items()]
    scores = pool.map(get_seq, ((n, s) for n, s in seq_list))
    pool.close()
    t2 = time.time()

    print("alignment finished in : " + str(round((t2-t1)/60, 2)))
    print("writing aligned file")
    if not os.path.isdir("alignedFile/"): os.mkdir("alignedFile/")

    alnRecord = []
    alnRecord_deleted = []
    for name, seq, ali_scr in scores:
        if len(seq) <= 184:
            alnRecord .append(SeqRecord(Seq(seq), id = name, name="", description=""))
            # nb_mut = sum(1 for s, az in zip(seq, ref_seq) if s != az and s != "-")
        else: 
            alnRecord_deleted.append(SeqRecord(Seq(seq), id = name, name="", description=""))
        # print(f">{name} SCORE={ali_scr} NB_MUT={nb_mut}\n{seq}")

        
    SeqIO.write(alnRecord, "alignedFile/"+ infile .replace(".fasta", "_aligned.fasta"), "fasta-2line")
    SeqIO.write(alnRecord, "alignedFile/"+ infile .replace(".fasta", "_aligned_deleted.fasta"), "fasta-2line")
    os.remove(infile)
    print("finish")