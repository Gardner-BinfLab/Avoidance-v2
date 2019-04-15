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
from itertools import cycle, islice
from functools import partial
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
    parser.add_argument('-s', '--sequence',
                        type=valid_file,
                        metavar='STR',
                        help='Sequences in fasta or csv format',
                        required='True')
    parser.add_argument('-U', '--utr5',
                        metavar='STR/INT',
                        help='Use an integer if 5UTR presence, e.g., -U 1. Use your own 5UTR sequence if your plasmid backbone is not of pET vector. Default = GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT')
    parser.add_argument('-x', '--execute',
                        help='Run RNAplfold multiprocessing',
                        action="store_true")
    parser.add_argument('-W', '--winsize',
                        default='210',
                        metavar='INT',
                        help='Average the pair probabilities over windows of given size. An RNAplfold option. Default = 210')                
    parser.add_argument('-u', '--ulength',
                        default='210',
                        metavar='INT',
                        help='Compute the mean probability that regions of length 1 to a given length are unpaired. An RNAplfold option. Default = 210') 
    parser.add_argument('-c', '--calculate',
                        help='Calculate sum of opening energy. Output .out, a tab delimited file',
                        action="store_true") 
    parser.add_argument('-i', '--ipos',
                        type=int,
                        metavar='INT',
                        default=71,
                        help='Position i of an input sequence. Related to option -c. e.g., if the length of 5UTR is 71 and l of choice is 43, to calculate sum of the opening energy of the first 30 nt of CDS use option -i 71 -j 101 -l 43. Default = 71')
    parser.add_argument('-j', '--jpos',
                        type=int,
                        default=101,
                        metavar='INT',
                        help='Position j of an input sequence. Related to option -c. Default = 101')
    parser.add_argument('-l', '--length',
                        type=int,
                        default=43,
                        metavar='INT',
                        help='l in _openen file. Related to option -c. Default = 43')
    parser.add_argument('-S', '--stack',
                        help='Stack _openen files to single columns and concatenate them as a pandas dataframe. Output .pkl, a pandas dataframe in pickle format',
                        action="store_true")
    parser.add_argument('-r', '--remove',
                        help='Remove _openen and .ps files',
                        action="store_true")
    parser.add_argument('-o', '--output',
                        metavar='STR',
                        default='openen', 
                        help='Output file name for .out and .pkl. Related to option -c and -S. Default = openen')
    parser.add_argument('-p', '--processes',
                        type=int,
                        metavar='INT',
                        help='Number of processes to spawn. Default = half of the number of CPU')

    results = parser.parse_args(args)
    return (results.sequence, results.execute, results.winsize, results.ulength, \
            results.calculate, results.ipos, results.jpos, results.length, \
            results.utr5, results.stack, results.remove, results.output, results.processes)


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
    timerthread = threading.Thread(target = time_count, args = ())
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
    fasta_df[['Accession','Sequence']]=fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['Accession'] = '>' + fasta_df['Accession']
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n','', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df.set_index('Accession',inplace=True)
    fasta_df = fasta_df[fasta_df.Sequence != '']
    final_df = fasta_df.dropna()
    return final_df



def openen(W, u, seq):
    w_par = '-W ' + str(W)
    u_par = '-u ' + str(u)
    run(['RNAplfold', w_par, u_par, '-O'], input=str.encode(seq))



def sum_openen(i, j, l, file):
    w = pd.read_csv(file, sep='\t', skiprows=2, header=None)[i:j][l].sum()
    return file.replace('_openen',''), w



def stack_openen(f):
    n = f.replace('_openen', '').split()
    n = pd.DataFrame(n)
    n.columns = ['id']#given that the 5UTR length is i
#sum opening energy from position i to j ntdefault=71, of CDS at l
#optimum at i=71, j=101, l=43
    d = pd.read_csv(f, sep='\t', skiprows=2, nrows=235, header=None)
    d = d.set_index(0).stack().to_frame()
    d = d[0].apply(lambda x: round(x, 4)).to_list()
    d = pd.DataFrame(d).T
    d = pd.concat([n,d], axis=1)
    return d
default=71,

        
def main():
    base,ext = os.path.splitext(s)
    if ext.lower() in ('.fasta', '.fa', '.fas', '.fna'):
        seq = fasta_to_dataframe(s).reset_index()
    else:
        seq = pd.read_csv(s, skiprows=1, header=None)

    startTime = datetime.now()

    label = (n for n in cycle(list(range(0,p))))
    label = pd.DataFrame({'label': list(islice(label, len(seq)))})
    df = pd.concat([label, seq], axis=1)
    _openen = (df['Accession'].str.replace('>', '') + '_openen').tolist()


    if x is True:
        print('\nRunning RNAplfold using', p, 'processes...')
        print_time()
        
        if U is None:
            utr5 = 'GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
            df['fasta'] = df['Accession'] + '\n' + utr5 + seq['Sequence'] + '\n'
            fasta = df.groupby('label')['fasta'].apply(''.join).tolist()
        else:
            try:
                int(u)
                df['fasta'] = df['Accession'] + '\n' + seq['Sequence'] + '\n'
                fasta = df.groupby('label')['fasta'].apply(''.join).tolist()
            except ValueError:     
                df['fasta'] = df['Accession'] + '\n' + u + seq['Sequence'] + '\n'
                fasta = df.groupby('label')['fasta'].apply(''.join).tolist()
        
        groups = len(fasta)
        p1 = Pool(p)
        results = []
        progress(0, groups)
        openen_func = partial(openen, W, u)
        for a in p1.imap_unordered(openen_func, fasta):
            results.append(a)
            progress(len(results), groups)
            
        _stop_timer.set()
        p1.close()
        p1.join()    
                
    else:
        print('\nSkipped RNAplfold (no option -x given)!', flush = True)  


    if c is True:
        print('\nCalculating sum of opening energy using', p, 'processes...')
        sum_func = partial(sum_openen, i, j, l)
        p2 = Pool(p)
        plfold =  p2.imap_unordered(sum_func, _openen)
        p2.close()
        p2.join()
        
        results = list(plfold)
        d = pd.DataFrame(results)
        d.columns = ['Accession', 'openen']    
        filename = o + '.out'
        d.to_csv(filename, sep='\t', index=False, encoding='utf-8')
    else:
        print('\nSkipped calculation (no option -c given)!', flush = True)  


    if S is True:
        print('\nParsing _openen using', p, 'processes...')
        print_time()

        d = pd.DataFrame()
        p3 = Pool(p)
        progress(0, len(_openen))
        for b in p3.imap_unordered(stack_openen, _openen):
            d = pd.concat([d, b], ignore_index=True)
            progress(len(d), len(_openen))
        
        _stop_timer.set()    
        p3.close()
        p3.join()
        
        filename = o + '.pkl'
        d.to_pickle(filename)
    
    else:
        print('\nSkipped parsing (no option -t given)!', flush = True)


    if r is True:
        print('\nRemoving temporary files...')
        ps = (df['Accession'].str.replace('>', '') + '_dp.ps').tolist()
        p4 = Pool(p)
        p4.imap_unordered(os.remove, _openen)
        p4.imap_unordered(os.remove, ps)
        p4.close()
        p4.join()
    
    else:
        print('\nSkipped removing _openen and .ps files (no option -r given)!', flush = True)


    print('\nTime taken', datetime.now() - startTime, flush = True)



if __name__ == "__main__":
    s, x, W, u, c, i, j, l, U, S, r, o, p = check_arg(sys.argv[1:])
    if U is None:
        utr5 = 'GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
        print('\n5UTR is', utr5, flush = True)
    else:
        try:
            int(U)
            print('\n5UTR presence', flush = True)
        except ValueError:
            print('\n5UTR is', U, flush = True)
                
    if p is None:
        p = int(multiprocessing.cpu_count()/2)
    main()
