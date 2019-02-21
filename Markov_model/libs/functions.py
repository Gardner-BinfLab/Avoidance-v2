import datetime
import sys
import RNA
import pandas as pd
import numpy as np
from libs import data
from numpy.random import choice


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


def fasta_to_dataframe(input_file):
    '''
    convert fasta to pandas dataframe
    '''
    sequence_df = pd.DataFrame(columns==[1,0])
    fasta = []
    test = []
    with open(input_file) as file:
        for line in file:
            line = line.strip()
            if not line:
               continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                sequence_df.loc[line,1] = line[1:]
                if active_sequence_name not in fasta:
                    test.append(''.join(fasta))
                    fasta = []
                continue
            sequence = line
            fasta.append(sequence)
    if fasta:
        test.append(''.join(fasta))
    for i, row in enumerate(test):
        sequence_df[0][i-1] = row
    return(sequence_df)




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
            print("Something's wrong with the input sequences.\n")
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


def score(sequence_df,prob_data,back_prob_data):
    '''
    bitscore a list of sequences using a foreground and a background model
    '''
    sequence_score=pd.DataFrame(columns=['scores'],index=list(range(len(sequence_df))))
    length_to_score = max([len(sequence_df[0][i]) for i in range(len(sequence_df))]) #score for full length
    for seq in range(len(sequence_df)):
        sequence = sequence_df[0][seq].lower()
        length = sequence_length(sequence) if len(sequence)<length_to_score\
                 else length_to_score
        codons = splitter(sequence,length)
        scores_df = pd.DataFrame(columns=['scores'],index=codons)
        back_scores_df = pd.DataFrame(columns=['scores'],index=codons)
        stop = ['tag','taa','tga']
        if bool(set(stop).intersection(codons[:len(codons)-1])) == False:
            scores_df = pd.DataFrame(columns=['scores'],index=codons)
            back_scores_df = pd.DataFrame(columns=['scores'],index=codons)
            for i in range(len(codons)):
                try:
                    if i < len(prob_data.columns)-1: #1 for average
                        scores_df.loc[codons[i],'scores'] = np.log2(prob_data.loc[codons[i],i])[0]
                        back_scores_df.loc[codons[i],'scores'] = np.log2(back_prob_data.loc[codons[i],i])[0]
                    else:
                        scores_df.loc[codons[i],'scores'] = np.log2(prob_data.loc[codons[i],'codon_prob'])[0]
                        back_scores_df.loc[codons[i],'scores'] = np.log2(back_prob_data.loc[codons[i],'codon_prob'])[0]


                except RuntimeWarning: #catch for log zero error
                    scores_df.loc[codons[i],'scores'] = np.log2(prob_data.loc[codons[i],'codon_prob'])[0]
                    back_scores_df.loc[codons[i],'scores'] = np.log2(back_prob_data.loc[codons[i],'codon_prob'])[0]
        else:
            print('\nStop codons encountered before the last position for sequence : ', seq)
            pass

        sequence_score.loc[seq,'scores'] = scores_df.sum()[0] - back_scores_df.sum()[0] 
        progress(seq,len(sequence_df))
        
    return sequence_score




def rna_ss(seq):
    '''calculates mean free energy of sequence using RNAlib
    '''
    ss,mfe = RNA.fold(seq)
    return mfe


def mutate(sequence,prob_df,length = 30):
    '''synonymously mutate sequence. The probability of codon from the model is used as a 
    weght to randomly pick the new codon.
    '''
    codons = functions.splitter(sequence,length)
    mutable_codon_position = choice(list(range(len(codons)))) #randomly choose pos to mutate
    mutable_synonymous_codons = data.aa2codon[data.codon2aa[codons[mutable_codon_position]]]
    
    #with probs of codon as weights for random picking
    try:
        probs_of_codons_for_emission = prob_df.loc[mutable_synonymous_codons,mutable_codon_position]
    except IndexError:
        probs_of_codons_for_emission = prob_df.loc[mutable_synonymous_codons,'codon_prob']
    
    probs_sum = sum(probs_of_codons_for_emission)
    mutated_codon = choice(mutable_synonymous_codons, p = [float(i)/probs_sum\
                                                        for i in probs_of_codons_for_emission])
    
    
    #without weights
    #mutated_codon = random.choice(mutable_synonymous_codons)
    new_seq = sequence[:mutable_codon_position*3]+ mutated_codon + sequence[mutable_codon_position*3+3:]
    return new_seq



def sim_anneal(sequence,prob_df,length = 30,niter=100):
    '''
    preforms a simulated annealing to maximize the secondary structure of sequence
    '''
    utr='ggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacat'
    seq = sequence[:length] 
    temp = np.linspace(1,0.001,niter)
    scurr = seq
    sbest = seq
    for i in range(niter):
        T = temp[i]
        snew = mutate(sbest,prob_df,length)
        if rna_ss(utr[-length:]+ snew) >= rna_ss(utr[-length:] + scurr):
                scurr = snew
                if rna_ss(utr[-length:]+scurr)>=rna_ss(utr[-length:]+sbest):
                    sbest = snew
        elif np.exp(-(rna_ss(utr[-length:]+scurr)-rna_ss(utr[-length:]+snew))/T) <= np.random.rand(1)[0]:
            scurr = snew
        functions.progress(i,niter)
    print('\nmfe is',rna_ss(utr[-length:]+sbest),sbest)
    annealed_seq = sbest + sequence[length:]
    return annealed_seq

