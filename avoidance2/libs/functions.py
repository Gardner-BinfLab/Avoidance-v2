import datetime
import sys
import pandas as pd
import numpy as np
from libs import data


def progress(iteration, total):   
    bars_string = int(float(iteration) / float(total) * 50.)
    sys.stdout.write(
        "\r[%-50s] %d%% (%s/%s)" % (
            '='*bars_string, float(iteration) / float(total) * 100,
            iteration,
            total
        ) 
    )
    sys.stdout.flush()
    if iteration == total:
        print('Finished!') 

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
 

def positionwise_codons(sequences,length_to_train,trainall=False):
    '''
    returns dictionary of positions(index) with codons.
    this will be converted to pandas dataframe soon
    ***deprecated. use codons_to_df**
    '''
    codon_dict={} 
    for sequence in sequences:
        length = sequence_length(sequence) if trainall==True or \
                 len(sequence)<length_to_train else length_to_train 
        codon_list = splitter(sequence.lower(),length)
        for i in range(len(codon_list)):
            try:
                codon_dict.setdefault(i, [])
                codon_dict[i].append(codon_list[i])
            except KeyError:
                codon_dict[i] = codon_list[i]

    return codon_dict

def codons_to_df(sequences,length_to_train,trainall=False):
    '''
    returns a dataframe of codons in given sequences
    position is column and sequence number is index
    '''
    codon_df = pd.DataFrame()
    index = 0
    skipped = 0
    for sequence in sequences:
        try:
            length = sequence_length(sequence) if trainall==True or \
                     len(sequence)<length_to_train else length_to_train
        except TypeError:
            print("something's wrong with the input sequences.\n")
            break
                  
        codon_list = splitter(sequence.lower().replace('u','t'),length)
        stop = ['tag','taa','tga']
        if bool(set(stop).intersection(codon_list[:len(codon_list)-1])) == False:
            for i in range(len(codon_list)):
                codon_df.at[index,i]=codon_list[i]
            
        else:
            print('\nStop codon encountered somewhere before the last position.')
            skipped += 1
        progress(index,len(sequences)-1)
        index += 1
    if skipped > 0:
        print(skipped, 'sequence(s) were skipped because we encountered stop codons.')
    

    return codon_df


def calc_cond_prob(codon_df):
    '''
    returns a first order conditional probability of codons
    '''
    cond_prob_df = []
    for i in range(len(codon_df.columns)-1): #conditional except last one
        cond_prob_df.append((codon_df.groupby([i, i+1]).count()\
                             / codon_df.groupby(i).count()))
        progress(i,len(codon_df.columns)-2)
    return cond_prob_df



    
def train(codon_df):
    '''
    uses the codon dataframe to calculate probabilities
    index is the codon name and column number is the position
    '''
    codons = [key for key,value in data.codon2aa.items()]
    prob_df = pd.DataFrame(index=[codons])
    for i in range(len(codon_df.columns)):
        prob_df[i]=np.nan*len(prob_df)
        probs=codon_df.groupby(i).size().\
                                      div(codon_df[i].count())
        for item in codon_df[i]:
            if str(item)!= 'nan':
                prob_df[i].loc[item]=probs.loc[item]
        progress(i,len(codon_df.columns)-1)    
    return prob_df






