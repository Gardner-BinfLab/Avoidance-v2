import os
import sys
import time
import pickle
import argparse
import numpy as np
import pandas as pd
from libs import functions





def check_arg(args=None):
    parser = argparse.ArgumentParser(description='This generates scores for given sequences.')
    parser.add_argument('-f', '--file',
                        help='Sequences in CSV. FASTA support coming soon!',
                        required='True')
    parser.add_argument('-m', '--model',
                        help='Trained model file in .pkl format.',
                        required='True')
    parser.add_argument('-b', '--background',
                        help='Background model file in .pkl format.',
                        required='True')
    parser.add_argument('-o', '--output',
                        help='Output file name',
                        default='scores')


    results = parser.parse_args(args)
    return (results.file,
            results.model,
            results.background,
            results.output)




def main():
    #os.makedirs(os.path.join(pathlib.Path.cwd(),'results','scores',''),exist_ok=True)
    mypath = os.path.join(os.getcwd(),'results','scores','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','scores',''))
    
    #mask NaNs with mean
    for i in range(prob_data.shape[1]-1):
        prob_data[i].fillna(prob_data['mean'], inplace=True)
    for i in range(back_prob_data.shape[1]-1):
        back_prob_data[i].fillna(back_prob_data['mean'], inplace=True)
    
    sequence_df = pd.read_csv(f,skiprows=1,header=None)
    sequence_score=pd.DataFrame(columns=['scores','scores_per_codon'],index=list(range(len(sequence_df))))
    length_to_score = max([len(sequence_df[0][i]) for i in range(len(sequence_df))]) #score for full length
    for seq in range(len(sequence_df)):
        sequence = sequence_df[0][seq].lower()
        length = functions.sequence_length(sequence) if len(sequence)<length_to_score\
                 else length_to_score
        codons = functions.splitter(sequence,length)
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
                        scores_df.loc[codons[i],'scores'] = np.log2(prob_data.loc[codons[i],'mean'])[0]
                        back_scores_df.loc[codons[i],'scores'] = np.log2(back_prob_data.loc[codons[i],'mean'])[0]


                except RuntimeWarning: #catch for log zero error
                    scores_df.loc[codons[i],'scores'] = np.log2(prob_data.loc[codons[i],'mean'])[0]
                    back_scores_df.loc[codons[i],'scores'] = np.log2(back_prob_data.loc[codons[i],'mean'])[0]
        else:
            print('\nStop codons encountered before the last position for sequence : ', seq)
            pass
        


        #sequence_score.loc[seq,'scores'] = scores_df.sum()[0] - back_scores_df.sum()[0]          
        sequence_score.loc[seq,'scores_per_codon'] = scores_df.mean()[0] - back_scores_df.mean()[0] 
        functions.progress(seq,len(sequence_df))


    print('\nExporting results..')    
    #sequence_df['scores'] = sequence_score['scores']
    sequence_df['scores'] = sequence_score['scores_per_codon']
    sequence_df.to_csv(mypath+'_'+o+'_'+time.strftime("%Y%m%d-%H%M%S") + '.csv'
                      ,sep=',', encoding='utf-8', index=False)


if __name__ == '__main__':
    f,m,b,o= check_arg(sys.argv[1:])
    with open(m, "rb") as model_file:
        prob_data = pickle.load(model_file)
    with open(b, "rb") as model_file:
        back_prob_data = pickle.load(model_file)
    main()

