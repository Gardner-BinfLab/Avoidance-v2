#!/usr/bin/env python

from builtins import range
import sys
import math
import argparse
from re import sub
import collections
import pandas as pd
from Bio import SeqIO
import RNA
import forgi.graph.bulge_graph as cgb


def sliding_windows(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break


def basepair(seqfile, win, step):
    for seq_record in SeqIO.parse(seqfile, "fasta"):
        seq = str(seq_record.seq)
        id = str(seq_record.id)
        for subseq in sliding_windows(seq, win, step):
            dotbracket = RNA.fold(subseq)[0]
            bg = cgb.BulgeGraph.from_dotbracket(dotbracket)
            bg.from_dotbracket(dotbracket)
            for s in bg.stem_iterator():
                for i in range(bg.stem_length(s)):
                    l=(bg.defines[s][0]+i)%3
                    r=(bg.defines[s][3]-i)%3
                    yield id, l, r


def count_basepair(seqfile, win, step):
    g = basepair(seqfile, win, step)
    for i in g:
        c = i[0]
        if (i[1]==0 and i[2]==2) or (i[1]==2 and i[2]==0) or (i[1]==1 and i[2]==1):
            yield c,'c1'
        elif (i[1]==0 and i[2]==1) or (i[1]==1 and i[2]==0) or (i[1]==2 and i[2]==2):
            yield c,'c2'
        else:
            yield c,'c3'



def main():
    counter = collections.Counter(count_basepair(args.i, args.w, args.s))
    df = pd.DataFrame.from_dict(counter, orient='index').reset_index()
    df = pd.concat([df['index'].apply(pd.Series), df[0]], axis = 1)
    df.columns = ['Accession', 'basepair', 'counts']
    df = df.pivot(index='Accession', columns='basepair', values='counts')
    df['Total']=df['c1']+df['c2']+df['c3']
    df['C1']=df['c1']/df['Total']*100
    df['C2']=df['c2']/df['Total']*100
    df['C3']=df['c3']/df['Total']*100
    df.to_csv(r'basepair.out',sep='\t')
    return


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description="Calculating RNA basepair classes")
    parser.add_argument('-i',type=str,metavar='STR',help="a FASTA file of CDS",required=True)
    parser.add_argument('-w',type=int,metavar='INT',help="Window size (nt)",required=True)
    parser.add_argument('-s',type=int,metavar='INT',help="Step (nt)",required=True)
    args=parser.parse_args()
    main()