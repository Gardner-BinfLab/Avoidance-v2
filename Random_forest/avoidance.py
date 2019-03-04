#!/usr/bin/env python

import os
import sys
import argparse
import itertools
import subprocess
import pandas as pd
import multiprocessing
from itertools import cycle
from libs import functions
from datetime import datetime
from multiprocessing import Pool
from subprocess import run, PIPE



def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta', '.fa', '.fas', '.fna'):
        raise argparse.ArgumentTypeError('File must have a fasta or csv extension')
    return param



def check_arg(args=None):
    parser = argparse.ArgumentParser(description='RNAup wrapper using multiprocesses')
    parser.add_argument('-m', '--mrna',
                        type=valid_file,
                        metavar='STR',
                        help='mRNA sequences in fasta or csv format',
                        required='True')
    parser.add_argument('-n', '--ncrna',
                        type=valid_file,
                        metavar='STR',
                        help='ncRNA sequences in fasta or csv format',
                        required='True')
    parser.add_argument('-l', '--length',
                        metavar='INT',
                        help='First N nt of mRNAs to calculate interactions for. Default = 30 nt')
    parser.add_argument('-o', '--output',
                        metavar='STR',
                        help='Output file name',
                        default='avoidance')
    parser.add_argument('-p', '--processes',
                        type=int,
                        metavar='INT',
                        help='Number of processes to spawn. Default = 16')

    results = parser.parse_args(args)
    return (results.mrna,
            results.ncrna,
            results.length,
            results.output,
            results.processes)



def interaction_calc(seq):
    proc = run(['RNAup', '-b', '-o', '--interaction_first'], stdout=PIPE,stderr=subprocess.DEVNULL,
               input=str.encode(seq)) #input is a stdin object so encode input str
    return str(proc.stdout).replace("\\n", " ").replace("b'", "")



def main():
    base,ext = os.path.splitext(m)
    if ext.lower() in ('.fasta', '.fa', '.fas', '.fna'):
        mrna = functions.fasta_to_dataframe(m)
    else:
        mrna = pd.read_csv(m, skiprows=1, header=None)
        
    base,ext = os.path.splitext(n)
    if ext.lower() in ('.fasta', '.fa', '.fas', '.fna'):
        ncrna = functions.fasta_to_dataframe(n)
    else:
        ncrna = pd.read_csv(n, skiprows=1, header=None)

    startTime = datetime.now()
    
    print('\nAssigning mRNA:ncRNA interactions...')
    startTime = datetime.now()
    mrna['mrna_seq'] = '>' + mrna[1].map(str) + '\n' + mrna[0].map(str).str[:length] + '\n'
    ncrna['ncrna_seq'] = '>' + ncrna[1].map(str) + '\n' + ncrna[0].map(str) + '\n'
    mrna_seq = [rows['mrna_seq'] for index,rows in mrna.iterrows()]   
    ncrna_seq = [rows['ncrna_seq'] for index,rows in ncrna.iterrows()]
    index = pd.MultiIndex.from_product([mrna_seq, ncrna_seq], names = ['mrna', 'ncrna'])
    sequence_df = pd.DataFrame(index=index).reset_index()
    df = sequence_df.pivot(index='mrna', columns='ncrna', values='ncrna')
    df['interaction_first'] = df.reset_index().values.sum(axis=1)
    print('\nWe took', datetime.now() - startTime, 'to assign these interactions!', flush=True)
  
    print('\nCalculating interactions using', p, 'processes...', flush=True)
    total_pairs = df.shape[0]
    my_pool = Pool(p)
    interactions = []
    functions.progress(0,total_pairs)
    for i in my_pool.imap_unordered(interaction_calc, df['interaction_first'], chunksize = int(total_pairs/p)):
        interactions.append(i)
        functions.progress(len(interactions), total_pairs)
        
    my_pool.close()
    my_pool.join()

    #parsing RNAup output
    seq_id = pd.Series(interactions).str.extractall(r'(>[\S]+)')[0].str.replace('>', '', regex=True).to_frame()
    mrna_id = (seq_id.loc[pd.IndexSlice[:, 0], :]).reset_index().set_index('level_0')
    ncrna_id = (seq_id.loc[pd.IndexSlice[:, 1:], :]).reset_index().set_index('level_0')
    binding_energy = pd.Series(interactions).str.extractall(r'(\(-[0-9]+\.[0-9]+)')[0].str.replace('(', '', regex=True).to_frame().reset_index().set_index('level_0')
    
    d = pd.concat([mrna_id.iloc[:,[1]], ncrna_id.iloc[:,[1]]], axis=1)
    match = cycle(list(range(len(ncrna))))
    d['match'] = [next(match) for i in range(len(d))]
    d = d.reset_index()
    d['level_0']=d['level_0'].astype(int)

    binding_energy = binding_energy.reset_index()
    binding_energy['level_0']=binding_energy['level_0'].astype(int)
    binding_energy['match']=binding_energy['match'].astype(int)
    
    d = pd.merge(binding_energy.reset_index(), d.reset_index(), on=['level_0', 'match'])
    d = d.iloc[:,[5,6,3]]
    d.columns = ['Accession', 'ncRNA', 'Binding_energy']

    filename = o + '.out'
    d.to_csv(filename, index=False, sep='\t', encoding='utf-8')

    print('\nWe took', datetime.now() - startTime, 'to finish the task!', flush = True)



if __name__ == "__main__":
    m,n,l,o,p = check_arg(sys.argv[1:])
    if p is None:
        p = 16
    if l is None:
        length = 30
        print('Calculations are of the first', length,'nt of mRNAs.', flush = True)
    else:
        try:
            length = int(l)
            print('Calculations are of the first', length,'nt of mRNAs.', flush = True)
        except ValueError:
            length = None
            print('Calculations are of full length.', flush = True)
    main()
