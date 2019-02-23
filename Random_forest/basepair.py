#!/usr/bin/env python

from __future__ import print_function
from builtins import range
import sys
import math
import argparse
from re import sub
from Bio import SeqIO
import RNA
import forgi.graph.bulge_graph as cgb


def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break


def basepair(seqfile):
    for seq_record in SeqIO.parse(seqfile, "fasta"):
        seq = str(seq_record.seq)
        id = str(seq_record.id)
        bg = cgb.BulgeGraph()
        for subseq in chunks(seq, 50, 5):
            dotbracket = RNA.fold(subseq)[0]
            bg.from_dotbracket(dotbracket)
            for s in bg.stem_iterator():
                for i in range(bg.stem_length(s)):
                    l=int((bg.defines[s][0]+i)%3)
                    r=int((bg.defines[s][3]-i)%3)
                    print(id, l, r, sep="\t")
    return


def main():
    basepair(args.i)
    return


if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="basepair.py")
    Argument_Parser.add_argument('-i',type=str,metavar='file',help="dotbracket with FASTA header.",required=True)
    args=Argument_Parser.parse_args()
    main()