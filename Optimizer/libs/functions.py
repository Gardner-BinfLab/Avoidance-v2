

import sys
import subprocess
import datetime
import pandas as pd
import numpy as np
from libs import data
from numpy.random import choice
from multiprocessing import Pool
from subprocess import run, PIPE 


def progress(iteration,total,message=None):
    if message is None:
        message = ''
    bars_string = int(float(iteration) / float(total) * 50.)
    print("\r|%-50s| %d%% (%s/%s) %s "% ('█'*bars_string+ "░" * \
                                     (50 - bars_string), float(iteration) / float(total) * 100,\
                                     iteration,total,message),end='\r',flush=True)

    if iteration == total:
        print('\nCompleted!') 


        
        
        
        
def read_fasta(f):
    fasta_df = pd.read_csv(f,sep='>', lineterminator='>',header=None)
    fasta_df[['accession','sequence']]=fasta_df[0].str.split('\n', 1, expand=True)
    fasta_df['accession'] = '>'+fasta_df['accession']
    fasta_df['sequence'] = fasta_df['sequence'].replace('\n','', regex=True)
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df.set_index('accession',inplace=True)
    fasta_df = fasta_df[fasta_df.sequence != '']
    final_df = fasta_df.dropna()
    return final_df
        
        
def sequence_length(seq):
    '''
    returns length of sequence in multiple of 3 by chopping extra positions
    '''
    if len(seq)%3 != 0:
        length = (len(seq)- len(seq)%3)
    else:
        length = len(seq)
    return length


def splitter(sequence,length):
    '''
    split sequence to codons
    '''
    length = length - length%3
    if length == 0:
        sys.stderr.write("Too small sequence. Trying for 3 nucleotides.\r")
        sys.stdout.flush()
        length = 3
    split_func = lambda sequence, n: [sequence[i:i+n] for\
                                    i in range(0, length, n)]
    return split_func(sequence,3)


def multiprocess_wrapper(function,input_dataframe):
    pools = Pool(input_dataframe.shape[0])
    pool_results = []
    for result in pools.imap_unordered(function, \
                                         input_dataframe):
        pool_results.append(result)



    pools.close()
    pools.join()
    return pool_results



def interaction_calc(seq):
    proc = run(['RNAup', '-b','-o'], stdout=PIPE,stderr=subprocess.DEVNULL,
               input=seq,encoding = 'utf-8')
    return str(proc.stdout).replace("\\n"," ").replace("b'","")




def rnaup_result_parser(raw_result_list,mrna_dataframe=None):
    '''return max interaction and a dataframe of all interactions
    '''
    #check if we are supplied a list (like for many sequences)
    #or single item like(for one seq vs several other seq)
    try:
        results_list = raw_result_list.copy()
    except AttributeError:
        results_list = [raw_result_list]
    
    
    interaction_df = pd.DataFrame({'unparsed_results':results_list})
    interaction_df[['accession','RNAup_output']] = interaction_df\
                                                   ['unparsed_results']\
                                                   .str.split(':break',1,\
                                                              expand=True)
    results_temp_df = pd.Series(interaction_df['RNAup_output']\
                                .str.extractall(r'((?<=\:).*?(?==))')[0])\
                                .str.split(pat='(', n=-1, expand=True)\
                                .drop(0, 1).astype(np.float64).unstack()
    if mrna_dataframe is not None:
        results_temp_df.columns = mrna_dataframe.index
    result_df = pd.concat([interaction_df, results_temp_df], axis=1)
    return results_temp_df,result_df



def rand_background(sequence,n=1000):
    '''random background
    '''
       
    length = sequence_length(sequence)
    codons = [k for k,v in data.codon2aa.items()]
    backgnd_seq = pd.DataFrame({'sequence':[''.join(np.random.choice(codons,length))\
                                            for _ in range(n)]})
          
    return backgnd_seq



def syn_background(sequence,n=1000):
    '''synonymous background
    '''
    sequence = sequence.lower()   
    length = sequence_length(sequence)
    codons = splitter(sequence,length)
    syn_seq = []
    for j in range(n):
        chain=''
        for i in range(len(codons)):
            possible_codons = data.aa2codon[data.codon2aa[codons[i]]]
            chain+=np.random.choice(possible_codons)
        syn_seq.append(chain)
    backgnd_seq = pd.DataFrame({'sequence':syn_seq})
          
    return backgnd_seq







def substitute_codon(sequence,num_of_subst=10):
    '''randomly substitute codons along the sequence at random positions
    '''
       
    length = sequence_length(sequence)
    #num_of_subst = np.random.choice(length)
    new_seq = sequence
    for i in range(num_of_subst):
        codons = splitter(new_seq,length)
        subst_codon_position = np.random.choice(list(range(len(codons)))) 
        subst_synonymous_codons = data.aa2codon[data.codon2aa[codons[subst_codon_position]]]
        subst_codon = np.random.choice(subst_synonymous_codons)
        new_seq = new_seq[:subst_codon_position*3]+ subst_codon +\
                    new_seq[subst_codon_position*3+3:]
          
    return new_seq


def hammingDistance(seq1,seq2):
    score_nt = 0

    for i in range(0,len(seq1)):
        if seq1[i] != seq2[i]:
            score_nt += 1
            
    return score_nt



