#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd
from libs import codon_usage, data


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.fasta', '.fa', '.fas', '.fna'):
        raise argparse.ArgumentTypeError('File must have a fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Calculating CAI and relative codon biosynthetic cost')
    parser.add_argument('-i', '--input',
                        type=valid_file,
                        metavar='STR',
                        help='CDS of interest in fasta format',
                        required='True')
    parser.add_argument('-r', '--reference',
                        type=valid_file,
                        metavar='STR',
                        help='reference CDS in fasta format',
                        required='True')
    parser.add_argument('-o', '--output',
                        metavar='STR',
                        help='Output file name',
                        default='scores')

    results = parser.parse_args(args)
    return (results.input,
            results.reference,
            results.output)


def fasta_to_dataframe(f):
    fasta_df = pd.read_csv(f,sep='>', lineterminator='>', header=None)
    fasta_df[['Accession','Sequence']] = fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['Accession'] = fasta_df['Accession']
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df.set_index('Accession', inplace=True)
    fasta_df = fasta_df[fasta_df.Sequence != '']
    final_df = fasta_df.dropna()
    return final_df


def splitter(sequence,length):
    '''
    split sequence to codons
    '''
    length = length - length%3
    if length == 0:
        sys.stderr.write("Too small sequence. Trying for 3 nucleotides.\r")
        sys.stdout.flush()
        length = 3
    split_func = lambda sequence, n: [sequence[i:i+n] for i in range(0, length, n)]
    return split_func(sequence,3)


def cost(seq):
    seq = seq.lower()
    given_seq = splitter(seq,len(seq))
    try:
        cost_values = [np.log(data.cost_table[codon]) for codon in given_seq[1:-1]] #except the start and stop codons
        score = np.exp(np.mean(cost_values))
    except KeyError:
        print('strange sequence or corrupted cost table!')
        return 0
    return score


def gc(seq):
    seq = seq.upper()
    given_seq = list(seq)
    d = pd.DataFrame(base for base in given_seq[3:-3])[0].value_counts().to_frame()
    score = d.loc[['G','C']][0].sum() / d[0].sum()
    return score


def gc3c(seq):
    seq = seq.upper()
    given_seq = splitter(seq,len(seq))
    try:
        d = pd.DataFrame(codon[2] for codon in given_seq[1:-1])[0].value_counts().to_frame()
        score = d.loc[['G','C']][0].sum() / d[0].sum()
    except KeyError:
        print('None of G or C are in the position 3 of all codons!')
        return 0
    return score


def run(seq, ref):
    heg = codon_usage.CodonAdaptationIndex()
    heg.generate_rscu(ref)
    for i in seq:
        yield gc(i), gc3c(i), codon_usage.CodonAdaptationIndex().cai_for_gene(i), heg.cai_for_gene(i), cost(i), i


def main():
    input = fasta_to_dataframe(i)

    d = pd.DataFrame(run(list(input['Sequence']), r)) # list(input['sequence'].map(str).str[3:-3])
    d.columns = ['GC','GC3C','CAI','CAI_HEG','Biosynthetic_cost','Sequence']
    d = pd.merge(input.reset_index(), d, on='Sequence').iloc[:,[0,2,3,4,5,6]]
    d = d.loc[d['Biosynthetic_cost'] != 0]
    filename = o + '.out'
    d.to_csv(filename, index=False, sep='\t', encoding='utf-8')



if __name__ == "__main__":
    i,r,o = check_arg(sys.argv[1:])
    main()


