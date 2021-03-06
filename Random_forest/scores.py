#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd
from libs import data, functions, features


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
#     parser.add_argument('-r', '--reference',
#                         type=valid_file,
#                         metavar='STR',
#                         help='reference CDS in fasta format',
#                         required='True')
    parser.add_argument('-o', '--output',
                        metavar='STR',
                        help='Output file name',
                        default='scores')

    results = parser.parse_args(args)
    return (results.input,
#             results.reference,
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


def cost_rna(seq):
    seq = seq.upper()
    given_seq = functions.splitter(seq,len(seq))
    excluded_codons = {'ATG', 'TGG', 'TGA', 'TAA', 'TAG'}
    codons = [codon for codon in given_seq if codon not in excluded_codons]
    seq = list(''.join(codons))
    try:
        cost_base = [data.cost_rna[base] for base in seq]
        score = np.sum(cost_base)
    except KeyError:
        print('strange sequence or corrupted ribonucleotide cost table!')
        return 0
    return score


def cost_protein(seq):
    seq = seq.upper()
    given_seq = functions.splitter(seq,len(seq))
    excluded_codons = {'ATG', 'TGG', 'TGA', 'TAA', 'TAG'}
    codons = [codon for codon in given_seq if codon not in excluded_codons]
    try:
        amino_acids = [data.codon2aa[codon] for codon in codons]
        cost = [data.cost_aa[aa] for aa in amino_acids]
        score = np.sum(cost)
    except KeyError:
        print('strange sequence or corrupted amino acid cost table!')
        return 0
    return score


def cost_codon(seq):
    seq = seq.upper()
    given_seq = splitter(seq,len(seq))
    excluded_codons = {'ATG', 'TGG', 'TGA', 'TAA', 'TAG'}
    codons = [codon for codon in given_seq if codon not in excluded_codons]
    try:
        cost_values = [np.log(data.cost_table[codon]) for codon in codons] #except the start and stop codons
        score = np.exp(np.mean(cost_values))
    except KeyError:
        print('strange sequence or corrupted cost table!')
        return 0
    return score


def gc(seq):
    seq = seq.upper()
    given_seq = functions.splitter(seq,len(seq))
    excluded_codons = {'ATG', 'TGG', 'TGA', 'TAA', 'TAG'}
    codons = [codon for codon in given_seq if codon not in excluded_codons]
    seq = list(''.join(codons))
    d = pd.DataFrame(base for base in seq)[0].value_counts().to_frame()
    score = d.loc[['G','C']][0].sum() / d[0].sum()
    return score


def gc3c(seq):
    seq = seq.upper()
    given_seq = splitter(seq,len(seq))
    excluded_codons = {'ATG', 'TGG', 'TGA', 'TAA', 'TAG'}
    codons = [codon for codon in given_seq if codon not in excluded_codons]
    try:
        d = pd.DataFrame(codon[2] for codon in codons)[0].value_counts().to_frame()
        score = d.loc[['G','C']][0].sum() / d[0].sum()
    except KeyError:
        print('None of G or C are in the position 3 of all codons!')
        return 0
    return score


def run(seq):
#     heg = codon_usage.CodonAdaptationIndex()
#     heg.generate_rscu(ref)
    for i in seq:
        yield gc(i), gc3c(i), features.Analyze(i).cai(), cost_codon(i), cost_rna(i), cost_protein(i), i


def main():
    input = fasta_to_dataframe(i)

    d = pd.DataFrame(run(list(input['Sequence'])))
    d.columns = ['GC','GC3C','CAI','Cost_codon','Cost_RNA','Cost_protein','Sequence']
    d = pd.merge(input.reset_index(), d, on='Sequence').iloc[:,[0,2,3,4,5,6,7]]
    d = d.loc[d['Cost_codon'] != 0]
    filename = o + '.out'
    d.to_csv(filename, index=False, sep='\t', encoding='utf-8')



if __name__ == "__main__":
    i,o = check_arg(sys.argv[1:])
    main()


