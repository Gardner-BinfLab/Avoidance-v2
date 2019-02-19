import sys
import os
import argparse
import pickle
import time
import pandas as pd
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
    
    mypath = os.path.join(os.getcwd(),'results','synonymous_sequences','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','synonymous_sequences',''))
    
    
    
    input_sequence = s.lower()
    syn_sequences = []
    for seq in range(n):
        length = functions.sequence_length(input_sequence)
        codons = functions.splitter(input_sequence,length)
        chain = ''



        for i in range(len(codons)):
            possible_codons = data.aa2codon[data.codon2aa[codons[i]]]
            try:
                probs_of_codons_for_emission = prob_data.loc[possible_codons,i]
                if probs_of_codons_for_emission.dropna().empty == True:
                    try:
                        probs_of_codons_for_emission = prob_data.loc[possible_codons,'mean']
                    except (IndexError,RuntimeWarning): #check if we are at the end
                        probs_of_codons_for_emission = prob_data.loc[possible_codons,'mean']
            except (KeyError,RuntimeWarning): #for positions longer then in the probability table
                probs_of_codons_for_emission = prob_data.loc[possible_codons,'mean']
                
                

            scores_sum = sum(probs_of_codons_for_emission)
            chain+=choice(possible_codons, p = [float(i)/scores_sum\
                                                for i in probs_of_codons_for_emission])

            #for random sequences
            #chain+=choice(prob_data.index,p=prob_data.loc[:,i])

        syn_sequences.append(chain.upper())
        functions.progress(seq,n)
    
    print("\nExporting sequences...")
    dataframe = pd.DataFrame(syn_sequences,columns=['synonymous_sequence'])
    filename = mypath+'_' + o +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    dataframe.to_csv(filename,sep=',', encoding='utf-8', index=False) 
        

        
if __name__ == '__main__':
    s,m,n,o= check_arg(sys.argv[1:])
    with open(m, "rb") as model_file:
        prob_data = pickle.load(model_file)
    main()

   
