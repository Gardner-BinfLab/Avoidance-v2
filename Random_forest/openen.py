#!/usr/bin/env python

import os
import sys
import glob
import time
import argparse
import itertools
import subprocess
import pandas as pd
import multiprocessing
import threading
from itertools import cycle, islice
from threading import Semaphore
from datetime import datetime
from multiprocessing import Pool
from subprocess import run, PIPE



def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta', '.fa', '.fas', '.fna'):
        raise argparse.ArgumentTypeError('File must have a fasta or csv extension')
    return param



def check_arg(args=None):
    parser = argparse.ArgumentParser(description='RNAplfold wrapper using multiprocesses')
    parser.add_argument('-i', '--input',
                        type=valid_file,
                        metavar='STR',
                        help='Sequences in fasta or csv format',
                        required='True')
    parser.add_argument('-u', '--utr5',
                        type=str,
                        metavar='STR',
                        default='GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT',
                        help='Use your own 5UTR sequence (71 nt) if your plasmid backbone is not of pET vector')
    parser.add_argument('-o', '--output',
                        metavar='STR',
                        default='openen', 
                        help='Output file name. Default = openen')
    parser.add_argument('-p', '--processes',
                        type=int,
                        metavar='INT',
                        help='Number of processes to spawn. Default = total number of CPU')

    results = parser.parse_args(args)
    return (results.input, results.utr5, results.output, results.processes)



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
    fasta_df = pd.read_csv(f,sep='>', lineterminator='>',header=None)
    fasta_df[['accession','sequence']]=fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['accession'] = '>'+fasta_df['accession']
    fasta_df['sequence'] = fasta_df['sequence'].replace('\n','', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df.set_index('accession',inplace=True)
    fasta_df = fasta_df[fasta_df.sequence != '']
    final_df = fasta_df.dropna()
    return final_df



def openen(seq):
    run(['RNAplfold', '-W 210', '-u 210', '-O'], input=str.encode(seq))



def openen43(files):
    for i in files:
        w = pd.read_csv(i, sep='\t', skiprows=2, header=None)[71:101][43].to_frame().apply(sum).values[0]
        yield i.replace('_openen',''), w


        
def main():
    base,ext = os.path.splitext(i)
    if ext.lower() in ('.fasta', '.fa', '.fas', '.fna'):
        seq = fasta_to_dataframe(i).reset_index()
    else:
        seq = pd.read_csv(i, skiprows=1, header=None)

    startTime = datetime.now()
    
    print('Calculating open energy using', p, 'processes...')
    print_time()
    
    label = (x for x in cycle(list(range(0,p))))
    label = pd.DataFrame({'label': list(islice(label, len(seq)))})
    df = pd.concat([label, seq], axis=1)
    df['fasta'] = df['accession'] + '\n' + u + seq['sequence'] + '\n'
    fasta = df.groupby('label')['fasta'].apply(''.join).tolist()
    
    groups = len(fasta)
    my_pool = Pool(p)
    results = []
    progress(0,groups)
    for j in my_pool.imap_unordered(openen, fasta):
        results.append(j)
        progress(len(results), groups)
        
    _stop_timer.set()
    my_pool.close()
    my_pool.join()
    
    #parsing RNAplfold output
    files = glob.glob('*_openen')
    ps = glob.glob('*.ps')
    d = pd.DataFrame(openen43(files))
    d.columns = ['Accession', 'openen43']

    filename = o + '.out'
    d.to_csv(filename, sep='\t', index=False, encoding='utf-8')
    
    for k in files:
        os.remove(k)
    for l in ps:
        os.remove(l)
        
    print('\nWe took', datetime.now() - startTime, 'to finish the task!', flush = True)



if __name__ == "__main__":
    i, u, o, p = check_arg(sys.argv[1:])
    if p is None:
        p = multiprocessing.cpu_count()
    main()