import sys
import os
import argparse
import time
import pickle
import pandas as pd
import numpy as np
from libs import functions


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Training script')
    parser.add_argument('-f', '--file',
                        type=valid_file,
                        help='Training sequences in CSV or FASTA.',
                        required='True')
    parser.add_argument('-o', '--output',
                        help='Output file name.',
                        default = 'model')
    parser.add_argument('-l','--length',
                        type=int,
                        help='we train upto this length in sequences. Default 96 nt',
                        default = 96)
    parser.add_argument('-a','--trainall',
                        help='pass to train for whole length. Overrides previous arg if passed both',
                        action="store_true")

    results = parser.parse_args(args)
    return (results.file,
            results.output,
            results.length,
            results.trainall)

def main():
    
    #fix seed
    np.random.seed(12345)
    
    mypath = os.path.join(os.getcwd(),'results','model','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','model',''))

    base,ext = os.path.splitext(f)
    if ext.lower in ('.fasta','.fa'):
        seq_df = functions.fasta_to_dataframe(f)
    else:
        seq_df = pd.read_csv(f,skiprows=1,header=None)
        
    
    print('Reading codons..')
    codon_df = functions.codons_to_df(seq_df[0],l,a)
    print('\nTraining started. It may take a while..')
    print('\nZeroth order')
    prob_df = functions.train(codon_df)
    
    #count number of codons
    codons, counts = np.unique(codon_df.values,return_counts=True)
    
    #add an extra column at the end
    prob_df['codon_prob'] = np.nan
    
    #fill with codon probs
    for i in range(len(codons)):
        prob_df['codon_prob'][codons[i]] = counts[i]/counts.sum()
        
        
   #mask NaNs with codon probs
    for i in range(prob_df.shape[1]-1):
        prob_df[i].fillna(prob_df['codon_prob'], inplace=True)
    
    #normalize
    for j in range(prob_df.shape[1]-1):
        prob_df[j] = prob_df[j]/prob_df[j].sum()
    
    
    #print('\nOne') #nfirst order is overkill. So not doing it currently
    #cond_prob_df = functions.calc_cond_prob(codon_df)
    #results = [prob_df,cond_prob_df]
    print('\nExporting model data to file..')
    filename = mypath+str(o)+'_'+time.strftime("%Y%m%d-%H%M%S")+'.pkl'
    with open(filename, "wb") as export_file:
        pickle.dump(prob_df, export_file)


if __name__ == '__main__':
    f,o,l,a= check_arg(sys.argv[1:])
    main()


