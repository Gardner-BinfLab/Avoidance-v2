#!/usr/bin/env python

import os
import sys
import argparse
import itertools
import subprocess
import pandas as pd
import multiprocessing
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



def progress(iteration, total):   
    bars_string = int(float(iteration) / float(total) * 50.)
    sys.stdout.write(
        "\r|%-50s| %d%% (%s/%s)" % (
            '█'*bars_string+ "░" * (50 - bars_string), float(iteration) / float(total) * 100,
            iteration,
            total
        ) 
    )
    sys.stdout.flush()
    if iteration == total:
        print(' Completed!') 


def fasta_to_dataframe(f):
    fasta_df = pd.read_csv(f,sep='>', lineterminator='>',header=None)
    fasta_df[['accession','sequence']]=fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['accession'] = '>'+fasta_df['accession']
    fasta_df['sequence'] = fasta_df['sequence'].replace('\n','', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df.set_index('accession',inplace=True)
    fasta_df = fasta_df[fasta_df.sequence != '']
    final_df = fasta_df.dropna()
    return final_df



def interaction_calc(seq):
    proc = run(['RNAup', '-b', '-o', '--interaction_first'], stdout=PIPE,stderr=subprocess.DEVNULL,
               input=str.encode(seq)) #input is a stdin object so encode input str
    return str(proc.stdout).replace("\\n", " ").replace("b'", "")



def main():
    base,ext = os.path.splitext(m)
    if ext.lower() in ('.fasta', '.fa', '.fas', '.fna'):
        mrna = fasta_to_dataframe(m)
    else:
        mrna = pd.read_csv(m, skiprows=1, header=None)
        
    base,ext = os.path.splitext(n)
    if ext.lower() in ('.fasta', '.fa', '.fas', '.fna'):
        ncrna = fasta_to_dataframe(n)
    else:
        ncrna = pd.read_csv(n, skiprows=1, header=None)

    startTime = datetime.now()
    
    print('\nAssigning mRNA:ncRNA interactions...')
    startTime = datetime.now()
    mrna['mrna_seq'] = mrna.index + '\n' + mrna['sequence'].map(str).str[:length] + '\n'
    ncrna['ncrna_seq'] = ncrna.index + ':break' + '\n' + ncrna['sequence'].map(str) + '\n'
    mrna_seq = [rows['mrna_seq'] for index,rows in mrna.iterrows()]   
    ncrna_seq = [rows['ncrna_seq'] for index,rows in ncrna.iterrows()]
    index = pd.MultiIndex.from_product([mrna_seq, ncrna_seq], names = ['mrna', 'ncrna'])
    sequence_df = pd.DataFrame(index=index).reset_index()
    df = sequence_df.pivot(index='mrna', columns='ncrna', values='ncrna')
    df['interaction_first'] = df.reset_index().values.sum(axis=1)
    print('\nWe took', datetime.now() - startTime, 'to assign these interactions!', flush=True)
  
    print('\nCalculating interactions using', p, 'processes...', flush=True)
    groups = df.shape[0]
    my_pool = Pool(p)
    interactions = []
    progress(0,groups)
    for i in my_pool.imap_unordered(interaction_calc, df['interaction_first'], chunksize=int(groups/p)):
        interactions.append(i)
        progress(len(interactions), groups)
        
    my_pool.close()
    my_pool.join()

    #parsing RNAup output
    ncrna_id = pd.Series(interactions).str.extractall(r'(>[\S]+:break)')[0].str.replace('[>:break]', '', regex=True).to_frame().reset_index()
    mrna_id = pd.Series(interactions).str.extractall(r'(>[\S]+)')[0].str.replace('[>]', '', regex=True).to_frame().loc[pd.IndexSlice[:, 0], :].reset_index()
    binding_energy = pd.Series(interactions).str.extractall(r'(\(-[0-9]+\.[0-9]+)')[0].str.replace('(', '', regex=True).to_frame().reset_index()
    d = pd.merge(ncrna_id, binding_energy, on=['level_0','match'])
    d = pd.merge(mrna_id, d, on='level_0').iloc[:, [2,4,5]]
    d.columns = ['Accession', 'ncRNA', 'binding_energy']

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
