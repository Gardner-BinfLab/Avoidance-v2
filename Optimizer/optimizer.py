import sys
import os
import argparse
import time
import pandas as pd
import numpy as np
from libs import functions,data,features


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Training script')
    parser.add_argument('-b', '--background',
                        type=valid_file,
                        help='background sequences in csv or fasta.',
                        required='True')
    parser.add_argument('-m', '--mrna',
                        type=valid_file,
                        help='mrna sequence to optimize',
                        required='True')
    parser.add_argument('-o', '--output',
                        help='Output file name.',
                        default = 'optimized_sequences')

    results = parser.parse_args(args)
    return (results.background,
            results.mrna,
            results.output)


def optimization(sequence):
    optimization = features.Optimize(sequence,cai_mean, cai_std,\
                                     gc_cont_mean,gc_cont_std,ss_mean,\
                                     ss_std, avd_mean, avd_std,1000)
    optimized_sequence = [optimization.simulated_anneal() \
                          for _ in range(2)]
    result_df = pd.DataFrame({sequence:optimized_sequence})

    return result_df



def main():
    np.random.seed(12345)
    
    mypath = os.path.join(os.getcwd(),'results','optimized_sequences','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','optimized_sequences',''))

    base,ext = os.path.splitext(b)
    if ext.lower() in ('.fasta','.fa'):
        backgnd_seq = functions.read_fasta(b)
    else:
        backgnd_seq = pd.read_csv(b)
        
    base,ext = os.path.splitext(m)
    if ext.lower() in ('.fasta','.fa'):
        mrna_df = functions.read_fasta(m)
    else:
        mrna_df = pd.read_csv(m)        
        
        
    #calculate features for the background sequences first
    print('calculating features for background sequences..', flush = True)
    functions.progress(0,4)
    backgnd_seq['features'] = backgnd_seq['sequence'].apply(lambda x:\
                                                            features.Analyze(x))
    print('cai done!')
    backgnd_seq['cai'] = backgnd_seq['features'].apply(lambda x : x.cai())
    functions.progress(1,4)
    print('gc_cont done!')
    backgnd_seq['gc_cont'] = backgnd_seq['features'].apply(lambda x : x.gc_cont())
    functions.progress(2,4)
    print('sec_str done!')
    backgnd_seq['sec_str'] = backgnd_seq['features'].apply(lambda x : x.sec_str())
    functions.progress(3,4)
    print('avoidance')
    backgnd_seq['avoidance'] = backgnd_seq['features'].apply(lambda x :\
                                                             x.avoidance())
    functions.progress(4,4)
    
    #for z scores
    cai_mean, cai_std = np.mean(backgnd_seq['cai'].values),\
                        np.std(backgnd_seq['cai'].values)
    gc_cont_mean, gc_cont_std = np.mean(backgnd_seq['gc_cont'].values),\
                                np.std(backgnd_seq['gc_cont'].values)
    ss_mean, ss_std = np.mean(backgnd_seq['sec_str'].values),\
                      np.std(backgnd_seq['sec_str'].values)
    avd_mean, avd_std = np.mean(backgnd_seq['avoidance'].values),\
                      np.std(backgnd_seq['avoidance'].values)
    
    
    #results = functions.multiprocess_wrapper(optimization,mrna_df['sequence'])
    #final = pd.concat([dfs for df in results], axis=1)
    
    
    count = 0
    opt_df = pd.DataFrame()
    for sequence in mrna_df['sequence'] :
        optimization = features.Optimize(sequence,cai_mean, cai_std,gc_cont_mean,\
                 gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,1000)
        sequence_count = mrna_df.loc[mrna_df['sequence']\
                                                       == sequence].index[0]
        print('optimizing for sequence :', sequence_count)
        optimized_sequence = [optimization.simulated_anneal() \
                              for _ in range(10)]
        opt_df[sequence_count]=optimized_sequence
        count+=1
        functions.progress(count,mrna_df.shape[0])
        
        
    filename = mypath + o +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    opt_df.to_csv(filename,sep=',', encoding='utf-8', index=False)

        
        

if __name__ == '__main__':
    b,m,o= check_arg(sys.argv[1:])
    main()
        
