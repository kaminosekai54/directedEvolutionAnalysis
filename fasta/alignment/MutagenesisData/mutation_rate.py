"""Compute 
"""
from align import read_fasta
from sys import argv


def main():
    infile = argv[1]
    ref_seq = "GTGGTAAAATCTGCCTAAACGGGGAAACTCTCACTGAGACAATCCCGTGCTAAATCAGCAATAGCTGTAAATGCCTAACGACTACACGGTAGACAACTCTAAGAGTTGAAGGTATAGTCTAAACTGCAAGGTGACTTGCAGATATCGG"
    len_seq = len(ref_seq)

    h_seq = read_fasta(infile)
    nb_seq = len(h_seq)
    seq_list = [list(s) for n, s in h_seq.items()]
    col_wise = list(map(list, zip(*seq_list)))

    for i, ni in enumerate(ref_seq):
        mut = sum(el != ni for el in col_wise[i])/nb_seq
        print(i+1, mut)


if __name__ == '__main__':
    main()
