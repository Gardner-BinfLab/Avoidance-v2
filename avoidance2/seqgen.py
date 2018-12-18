import argparse
import sys
import os
import time
import pandas as pd
import progressbar
import numpy as np
from hmmlearn import hmm
from sklearn.externals import joblib
from itertools import zip_longest
from libs import data
from libs import functions



def check_arg(args=None):
    parser = argparse.ArgumentParser(description='This attempts to generate synonymous sequences.')
    parser.add_argument('-s', '--sequence',
                        help='Sequence to optimize',
                        required='True')
    parser.add_argument('-m', '--model',
                        help='HMM file in .pkl format.',
                        required='True')
    parser.add_argument('-n',type=int, 
                        help='Number of synonymous sequences to generate.',
                        required='True')
    parser.add_argument('-v',
                        help='Verbosity. Leave blank for no verbosity.',
                        action="store_true")

    results = parser.parse_args(args)
    return (results.sequence,
            results.model,
            results.n,
            results.v)



def main():
    #os.makedirs(os.path.join(pathlib.Path.cwd(),'results','syn_seq',''),exist_ok=True)
    mypath = os.path.join(os.getcwd(),'results','syn_seq','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','syn_seq',''))


    base_model = joblib.load(m)
    pbar = progressbar.NullBar(min_value=0, max_value=None)
    if v == True:
        print('\nGenerating sequences..')
        pbar = progressbar.ProgressBar(widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ',
                                            progressbar.ETA()], maxval=n) 
       

    scores=[]
    synonymous_sequence=[]
    pbar.start()
    for seq in range(n):

        input_sequence = s.lower()
        split_sequence = lambda x, n: [x[i:i+n] for i in range(0, len(x), n)]
        codons = split_sequence(input_sequence,3)
        new_path = base_model.sample(len(codons))[1]


        choosen_codon_numeric=[]
        for i in range(len(codons)):
            possible_codons = data.aa2codon[data.codon2aa[codons[i]]]
            numeric_codons = [data.codon_to_n[item] for item in possible_codons]
            scores_of_codons_for_emission=[base_model.emissionprob_[new_path[i]][item]
                                           for item in numeric_codons]
            #check for all zero values
            if all(scores == 0 for scores in scores_of_codons_for_emission)== True:
                raise ValueError("""Insufficient data to determine codon at this position. Using \
                data from previous position.""")
                scores_of_codons_for_emission=[base_model.emissionprob_[new_path[i]-1][item]
                                           for item in numeric_codons]
                
            choosen_codon_numeric.append(numeric_codons[
                scores_of_codons_for_emission.index(max(scores_of_codons_for_emission))])

        translated_codons = [data.n_to_codon[item] for item in choosen_codon_numeric]
        synonymous_sequence.append((''.join(translated_codons)).upper())
        scores.append(base_model.score(np.array(choosen_codon_numeric).reshape(-1,1)))
        
        #print(synonymous_sequence,base_model.score(np.array(choosen_codon_numeric).reshape(-1,1)))

        pbar.update(seq)
    pbar.finish()
    progressbar.streams.flush()  
    
        
    seq_with_scores = pd.DataFrame(list(zip_longest(synonymous_sequence,scores )),\
                          columns=['synonymous_sequence','scores'])
    print('Exporting scores to csv..')
    seq_with_scores.to_csv(mypath+'synonymous_sequence'+time.strftime("%Y%m%d-%H%M%S") + '.csv'
                      ,sep=',', encoding='utf-8', index=False)        
    

if __name__ == '__main__':
    s,m,n,v= check_arg(sys.argv[1:])
    main()
















