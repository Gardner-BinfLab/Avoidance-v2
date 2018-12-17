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
from itertools import zip_longest
#from threading import Thread


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='This generates scores for given sequences.')
    parser.add_argument('-f', '--file',
                        help='Sequences in CSV. FASTA support coming soon!',
                        required='True')
    parser.add_argument('-m', '--model',
                        help='HMM file in .pkl format.',
                        required='True')
    parser.add_argument('-v',
                        help='Verbosity. Leave blank for no verbosity.',
                        action="store_true")

    results = parser.parse_args(args)
    return (results.file,
            results.model,
            results.v)




def main():
    os.makedirs(os.path.join(pathlib.Path.cwd(),'results','scores',''),exist_ok=True)
    mypath =  os.path.join(pathlib.Path.cwd(),'results','scores','')


    df = pd.read_csv(f,skiprows=1,header=None)
    codon_list = functions.codon_splitter(df[0],3,v)
    lengths = [len(codon_list[index]) for index in range(len(codon_list))]
    codon_to_num = functions.create_mapping(codon_list,'encode')
    num_to_codon = functions.create_mapping(codon_list,'decode')
    X = functions.convert_to_hmm_data(codon_list,codon_to_num,v)
    base_model = joblib.load(m)
    scores = functions.calc_scores(base_model,X,lengths,v)


    scores_all = pd.DataFrame(list(zip_longest(df[0],scores )),\
                          columns=['sequences','scores'])
    #functions.helicopter('stop')
    print('Exporting scores to csv..')
    scores_all.to_csv(mypath+'scores_'+time.strftime("%Y%m%d-%H%M%S") + '.csv'
                      ,sep=',', encoding='utf-8', index=False)



if __name__ == '__main__':
    f,m,v= check_arg(sys.argv[1:])
    #if v == False:
    #    Thread(target=functions.helicopter('start')).start
    #    Thread(target=main()).start

    main()

