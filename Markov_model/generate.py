import sys
import os
import argparse
import pickle
import time
import pandas as pd
import numpy as np
from numpy.random import choice
from libs import functions
from libs import data


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='This attempts to generate synonymous sequences.')
    parser.add_argument('-s', '--sequence',
                        help='Sequence to optimize',
                        required='True')
    parser.add_argument('-m', '--model',
                        help='model file in .pkl format.',
                        required='True')
    parser.add_argument('-n',type=int, 
                        help='Number of synonymous sequences to generate.',
                        default=1000)
    parser.add_argument('-o', '--output',
                        help='Output file name.',
                        default = 'synonymous_seq')

    
    results = parser.parse_args(args)
    return (results.sequence,
            results.model,
            results.n,
            results.output)

    
def main():
    
    #fix seed
    np.random.seed(12345)
    
    mypath = os.path.join(os.getcwd(),'results','synonymous_sequences','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','synonymous_sequences',''))
    
    
    
    sequence = s.lower()
    #remove start and end
    input_sequence = sequence[3:-3]
    poss_term = 'ttttt'
    syn_sequences = []
    while True:
        length = functions.sequence_length(input_sequence)
        codons = functions.splitter(input_sequence,length)
        chain = ''



        for i in range(len(codons)):
            possible_codons = data.aa2codon[data.codon2aa[codons[i]]]
            try:
                probs_of_codons_for_emission = prob_data.loc[possible_codons,i]
                if probs_of_codons_for_emission.dropna().empty == True:
                    try:
                        probs_of_codons_for_emission = prob_data.loc[possible_codons,'codon_prob']
                    except (IndexError,RuntimeWarning): #check if we are at the end
                        probs_of_codons_for_emission = prob_data.loc[possible_codons,'codon_prob']
            except (KeyError,RuntimeWarning): #for positions longer then in the probability table
                probs_of_codons_for_emission = prob_data.loc[possible_codons,'codon_prob']
                
                

            scores_sum = sum(probs_of_codons_for_emission)
            chain+=choice(possible_codons, p = [float(i)/scores_sum\
                                                for i in probs_of_codons_for_emission])

            #for random sequences
            #chain+=choice(prob_data.index,p=prob_data.loc[:,i])
        final_seq = sequence[:3] + chain + sequence[-3:]
        if 'ttttt' not in final_seq and\
        'cacctgc' not in final_seq and \
        'gcaggtg' not in final_seq and \
        'ggtctc' not in final_seq and\
        'gagacc' not in final_seq and\
        'cgtctc' not in final_seq and\
        'gagacg' not in final_seq:
            syn_sequences.append(final_seq.upper())
            functions.progress(len(syn_sequences),n)

        
        
        if len(syn_sequences) == n:
            break
    
    print("\nExporting sequences...")
    dataframe = pd.DataFrame(syn_sequences,columns=['sequence'])
    filename = mypath + o +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    dataframe.to_csv(filename,sep=',', encoding='utf-8', index=False) 
        

        
if __name__ == '__main__':
    s,m,n,o= check_arg(sys.argv[1:])
    with open(m, "rb") as model_file:
        prob_data = pickle.load(model_file)
    main()

   
