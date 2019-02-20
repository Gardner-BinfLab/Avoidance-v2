import os
import sys
import time
import pickle
import argparse
import numpy as np
import pandas as pd
from libs import functions

def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param




def check_arg(args=None):
    parser = argparse.ArgumentParser(description='This generates scores for given sequences.')
    parser.add_argument('-f', '--file',
                        help='Sequences in CSV or FASTA',
                        type=valid_file,
                        required='True')
    parser.add_argument('-m', '--model',
                        help='Trained model file in .pkl format.',
                        required='True')
    parser.add_argument('-b', '--background',
                        help='Background model file in .pkl format. Default is uniform background.')
    parser.add_argument('-o', '--output',
                        help='Output file name',
                        default='scores')


    results = parser.parse_args(args)
    return (results.file,
            results.model,
            results.background,
            results.output)




def main():
    
    #fix seed
    np.random.seed(12345)
    
    
    #os.makedirs(os.path.join(pathlib.Path.cwd(),'results','scores',''),exist_ok=True)
    mypath = os.path.join(os.getcwd(),'results','scores','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','scores',''))

   #mask NaNs with mean
    #for i in range(prob_data.shape[1]-1):
    #    prob_data[i].fillna(prob_data['mean'], inplace=True)
    #for i in range(back_prob_data.shape[1]-1):
    #    back_prob_data[i].fillna(back_prob_data['mean'], inplace=True)


    base,ext = os.path.splitext(f)
    if ext.lower() in ['.fasta','.fa']:
        sequence_df = functions.fasta_to_dataframe(f)
    else:
        sequence_df = pd.read_csv(f,skiprows=1,header=None)

        
        
    sequence_score = functions.score(sequence_df,prob_data,back_prob_data)
    print('\nExporting results..')
    sequence_df.index = sequence_score.index
    sequence_df['scores'] = sequence_score['scores']
    sequence_df.to_csv(mypath+o+'_'+time.strftime("%Y%m%d-%H%M%S") + '.csv'
                      ,sep=',', encoding='utf-8', index=False)


if __name__ == '__main__':
    f,m,b,o= check_arg(sys.argv[1:])
    with open(m, "rb") as model_file:
        prob_data = pickle.load(model_file)
    if b is None:
        print('Using uniform background..')
        back_prob_data = pd.DataFrame(1/prob_data.shape[0], index=prob_data.axes[0],\
                                      columns=prob_data.columns)
    else:
        with open(b, "rb") as model_file:
            back_prob_data = pickle.load(model_file)
    main()

