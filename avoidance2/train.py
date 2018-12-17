import argparse
import sys
import os
import pathlib
import pandas as pd
import numpy as np
import time
from sklearn.externals import joblib
from hmmlearn import hmm
from libs import functions

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Training script')
    parser.add_argument('-f', '--file',
                        help='Training sequences in CSV. FASTA support coming soon!',
                        required='True')
    parser.add_argument('-n',type=int,
                        help='The number of hidden states to use. Defaults to 15',
                        default='15',
                        required='True')
    parser.add_argument('-i',type=int,
                        help='The number of iterations to perform for EM. Defaults to 500',
                        default='500',
                        required='True')
    parser.add_argument('-v',
                        help='Verbosity. Leave blank for no verbosity.',
                        action="store_true")

    results = parser.parse_args(args)
    return (results.file,
            results.n,
            results.i,
            results.v)


def main():
    if v == True:
        print('Reading input sequence file..')


    os.makedirs(os.path.join(pathlib.Path.cwd(),'results','hmm',''),exist_ok=True)
    mypath = os.path.join(pathlib.Path.cwd(),'results','hmm','')


    df = pd.read_csv(f,skiprows=1,header=None)
    codon_list = functions.codon_splitter(df[0],3,v)
    lengths = [len(codon_list[index]) for index in range(len(codon_list))]
    codon_to_num = functions.create_mapping(codon_list,'encode')
    num_to_codon = functions.create_mapping(codon_list,'decode')
    X = functions.convert_to_hmm_data(codon_list,codon_to_num,v)
    base_model = hmm.MultinomialHMM(n_components=n, n_iter=i, verbose = v)
    print('Training begins. This may take a while..')
    base_model.fit(X,lengths)
    print('Training completed.\nDumping model.')
    joblib.dump(base_model,mypath+'hmm_'+time.strftime("%Y%m%d-%H%M%S") + '.pkl')





if __name__ == '__main__':
    f,n,i,v= check_arg(sys.argv[1:])
    main()


