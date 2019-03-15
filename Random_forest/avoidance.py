#!/usr/bin/env python

import os
import sys
import time
import argparse
import itertools
import subprocess
import pandas as pd
import multiprocessing
import threading
from threading import Semaphore
from itertools import cycle, islice
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
    parser.add_argument('-b', '--batch',
                        help='Switch on batch processing for one ncRNA against many mRNAs',
                        action="store_true")
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
            results.batch,
            results.length,
            results.output,
            results.processes)


screen_lock = Semaphore(value=1)
_stop_timer = threading.Event() #global var for thread
def time_count():
    starttime = datetime.now()
    while not _stop_timer.isSet():
        screen_lock.acquire()
        time_message = '\r' + str(datetime.now() - starttime)
        sys.stdout.write(time_message)
        sys.stdout.flush()
        screen_lock.release()
        time.sleep(0.01)


def print_time():
    timerthread = threading.Thread(target = time_count,args = ())
    timerthread.start()

    
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
    fasta_df = pd.read_csv(f,sep='>', lineterminator='>', header=None)
    fasta_df[['accession','sequence']] = fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['accession'] = '>' + fasta_df['accession']
    fasta_df['sequence'] = fasta_df['sequence'].replace('\n', '', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df.set_index('accession', inplace=True)
    fasta_df = fasta_df[fasta_df.sequence != '']
    final_df = fasta_df.dropna()
    return final_df



def interaction_calc(seq):
    proc = run(['RNAup', '-b', '-o', '--interaction_first'], stdout=PIPE, stderr=subprocess.DEVNULL,
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
    
    if b is True: #one ncRNA vs multi mRNAs
        ncrna = ncrna.reset_index().reset_index().rename(index=str, columns={'index': 'interaction'})
               
        label = (x for x in cycle(list(range(0,p))))
        label = pd.DataFrame({'label': list(islice(label, len(mrna)))})
        df = pd.concat([label, mrna.reset_index()], axis=1)
        df['mrna'] = df['accession'] + ':break' + '\n' + df['sequence'].map(str).str[:length] + '\n'
        df = df.groupby('label')['mrna'].apply(''.join).to_frame().reset_index()
        df.insert(2, 'interaction', 0)
        df = pd.merge(df, ncrna, on='interaction')
        df['fasta'] = df['accession'] + '\n' + df['sequence'] + '\n' + df['mrna']
        fasta = df['fasta'].tolist()
        print('\nWe took', datetime.now() - startTime, 'to assign these interactions!', flush=True)

        print('\nCalculating interactions using', p, 'processes...', flush=True)
        print_time()
        my_pool = Pool(p)
        interactions = []
        progress(0, len(fasta))
        for i in my_pool.imap_unordered(interaction_calc, fasta):
            interactions.append(i)
            progress(len(interactions), len(fasta))

        _stop_timer.set()
        my_pool.close()
        my_pool.join()        
        

        
    else: #multiple ncRNAs and mRNAs
        ncrna['ncrna_seq'] = ncrna.index.map(str) + '\n' + ncrna['sequence'].map(str) + '\n'
        mrna['mrna_seq'] = mrna.index.map(str) + ':break' + '\n' + mrna['sequence'].map(str).str[:length] + '\n'
        mrna_seq = [rows['mrna_seq'] for index,rows in mrna.iterrows()]
        ncrna_seq = [rows['ncrna_seq'] for index,rows in ncrna.iterrows()]
        index = pd.MultiIndex.from_product([ncrna_seq, mrna_seq], names = ['ncrna', 'mrna'])
        sequence_df = pd.DataFrame(index = index).reset_index()
        df = sequence_df.pivot(index='ncrna',columns='mrna',values='mrna')
        df['interaction_first'] = df.reset_index().values.sum(axis=1)
        print('\nWe took', datetime.now() - startTime, 'to assign these interactions!', flush=True)

        print('\nCalculating interactions using', p, 'processes...', flush=True)
        print_time()
        groups = df.shape[0]
        my_pool = Pool(p)
        interactions = []
        progress(0,groups)
        for i in my_pool.imap_unordered(interaction_calc, df['interaction_first']):
            interactions.append(i)
            progress(len(interactions), groups)

        _stop_timer.set()
        my_pool.close()
        my_pool.join()

    #parsing RNAup output
    mrna_id = pd.Series(interactions).str.extractall(r'(>[\S]+:break)')[0].str.replace('>', '').str.replace(':break', '').to_frame().reset_index()
    ncrna_id = pd.Series(interactions).str.extractall(r'(>[\S]+)')[0].str.replace('>', '').to_frame().loc[pd.IndexSlice[:, 0], :].reset_index()
    binding_energy = pd.Series(interactions).str.extractall(r'(\([\S]+\.[0-9]+)')[0].str.replace('(', '').to_frame().reset_index()
    d = pd.merge(mrna_id, binding_energy, on=['level_0','match'])
    d = pd.merge(ncrna_id, d, on='level_0').iloc[:, [2,4,5]]
    d.columns = ['ncRNA', 'Accession', 'binding_energy']
    d = d.pivot(index='Accession',columns='ncRNA',values='binding_energy')
    filename = o + '.out'
    d.to_csv(filename, sep='\t', encoding='utf-8')

    print('\nWe took', datetime.now() - startTime, 'to finish the task!', flush = True)



if __name__ == "__main__":
    m,n,b,l,o,p = check_arg(sys.argv[1:])
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
            print('Calculations are of full-length.', flush = True)
    main()
