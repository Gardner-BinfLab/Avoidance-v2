import os
import sys
import re
import argparse
import itertools
import subprocess
import pandas as pd
import multiprocessing
from collections import deque
from datetime import datetime
from multiprocessing import Pool
from subprocess import Popen, PIPE


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Avoidance calculator script')
    parser.add_argument('-m', '--mrna',
                        type=valid_file,
                        help='mrna in fasta format',
                        required='True')
    parser.add_argument('-n', '--ncrna',
                        type=valid_file,
                        help='ncrna in fasta format',
                        required='True')
    parser.add_argument('-l','--length',
                        help='length to calculate interactions for. default = 30 nt')
    parser.add_argument('-o', '--output',
                        help='Output file name.',
                        default = 'avoidance')
    parser.add_argument('-p','--processes',
                        type=int,
                        help='number of process to spawn. Default = num cores/2')


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



def interaction_calc(seq_df):
    proc = Popen(['RNAup', '-b','-o','--interaction_first'],\
                 stdin = PIPE, stdout=PIPE, stderr=subprocess.DEVNULL)
    results = [str(proc.communicate(seq)[0]).replace("\\n"," ").replace("b'","") \
               for seq in pd.Series(seq_df)]

    return results



def main():
    
    
    mypath = os.path.join(os.getcwd(),'results','avoidance','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','avoidance',''))
        
        
    startTime = datetime.now()
    
    

    mrna = fasta_to_dataframe(m)
    ncrna = fasta_to_dataframe(n)                    
    
    
    temp_df = pd.DataFrame()
    temp_df['input_mrna'] = ncrna.index+'\n' + ncrna['sequence']
    ncrna_input = '\n'.join(temp_df['input_mrna'] )
    
    mrna['input'] = mrna.index+ ':break' + '\n'+ mrna['sequence'][:length] + '\n'+ ncrna_input
    mrna['input_encoded'] = mrna['input'].apply(lambda x: str.encode(x))
    
    
    total_seq = mrna.shape[0]
    chunks = 1 if p >= total_seq else p


    
    
    print("\nspawning ", p, " processes..", flush=True)
    my_pool = Pool(p)
    pool_results = deque()
    
    print("\ncalculating interactions.. this may take a while..")
    print('progressbar is updated after completion of every ',
    progress(0,total_seq)
    for result in my_pool.imap_unordered(interaction_calc, \
                                         mrna['input_encoded'],\
                                         chunksize = p):
        pool_results.append(result)
        completed = len(list(itertools.chain.from_iterable(pool_results)))
        progress(completed,total_seq) 
        
    my_pool.close()
    my_pool.join()
    
    print("we took ", datetime.now() - startTime, " to finish the task!",\
          flush = True)
    print('exporting results...', flush=True)
    
    
    interactions = list(itertools.chain.from_iterable(pool_results))
    interaction_df = pd.DataFrame({'results':interactions})
    interaction_df[['accession','RNAup_result']] = interaction_df['results'].str.split\
                                               (':break',1,expand=True)
    
    interaction_df['interaction'] = interaction_df['RNAup_result'].str.\
                                    extractall(r'(\(-[0-9]+\.[0-9]+)').unstack().\
                                    apply(','.join, 1).apply(lambda x: x.replace('(',''))
    
    filename = mypath + o +'_'+str(datetime.now()).replace(" ","_")+'.tsv'
    interaction_df.to_csv(filename, columns = ['accession','RNAup_result','interaction'],sep='\t',\
                          encoding='utf-8',index=None)
    
    print('done!')
    













if __name__ == "__main__":
    m,n,l,o,p = check_arg(sys.argv[1:])
    if p is None:
        #p = multiprocessing.cpu_count()
        p = int(multiprocessing.cpu_count()/2)
    if l is None:
        length = 30
        print('\ncalculations are for first ', length,' nucleotides.', flush = True)
    else:
        try:
            length = int(l)
            print('\ncalculations are for first ', length,' nucleotides.', flush = True)
        except ValueError:
            length = None
            print('\ncalculations are for full length.', flush = True)
    main()
