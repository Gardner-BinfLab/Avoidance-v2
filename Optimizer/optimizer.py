import sys
import os
import argparse
import time
import pandas as pd
import numpy as np
from libs import functions,data,features
from functools import partial


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Training script')
    parser.add_argument('-m', '--mrna',
                        type=valid_file,
                        help='mrna sequence to optimize in csv or fasta.',
                        required='True')
    parser.add_argument('-o', '--output',
                        help='Output file name.',
                        default = 'optimized_sequences')

    results = parser.parse_args(args)
    return (results.mrna,
            results.output)


def optimization(sequence):
    optimization = features.Optimize(sequence,cai_mean, cai_std,\
                                     gc_cont_mean,gc_cont_std,ss_mean,\
                                     ss_std, avd_mean, avd_std,1000)
    optimized_sequence = [optimization.simulated_anneal() \
                          for _ in range(2)]
    result_df = pd.DataFrame({sequence:optimized_sequence})

    return result_df


def dummy_func(optimization_df):
    optimization_df['annealed'] = optimization_df['objects'].apply(lambda x:x.simulated_anneal())
    return optimization_df



def main():
    np.random.seed(12345)
    
    mypath = os.path.join(os.getcwd(),'results','optimized_sequences','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','optimized_sequences',''))


    base,ext = os.path.splitext(m)
    if ext.lower() in ('.fasta','.fa'):
        mrna_df = functions.read_fasta(m)
    else:
        mrna_df = pd.read_csv(m)        
        
    backgnd_seq = functions.background(mrna_df['sequence'][0],10) 
    
    #calculate features for the background sequences first
    print('1) calculating features for background sequences..', flush = True)
    backgnd_seq['features'] = backgnd_seq['sequence'].apply(lambda x:\
                                                            features.Analyze(x))

    print('cai        ',end='\r',flush=True)
    backgnd_seq['cai'] = backgnd_seq['features'].apply(lambda x : x.cai())
    print('cai done!   ',end='\r',flush=True)
    
    print('gc_cont      ',end='\r',flush=True)
    backgnd_seq['gc_cont'] = backgnd_seq['features'].apply(lambda x : x.gc_cont())
    print('gc_cont done!',end='\r',flush=True)
    
    print('sec_str      ',end='\r',flush=True)
    backgnd_seq['sec_str'] = backgnd_seq['features'].apply(lambda x : x.sec_str())
    print('sec_str done!',end='\r',flush=True)
    
    print('avoidance     ',end='\r',flush=True)
    backgnd_seq['avoidance'] = backgnd_seq['features'].apply(lambda x :\
                                                             x.avoidance())
    print('avoidance done!',end='\r',flush=True)
    print('2) running optimizations. this may take a while...')
    
    
    #for z scores
    cai_mean, cai_std = np.mean(backgnd_seq['cai'].values),\
                        np.std(backgnd_seq['cai'].values)
    gc_cont_mean, gc_cont_std = np.mean(backgnd_seq['gc_cont'].values),\
                                np.std(backgnd_seq['gc_cont'].values)
    ss_mean, ss_std = np.mean(backgnd_seq['sec_str'].values),\
                      np.std(backgnd_seq['sec_str'].values)
    avd_mean, avd_std = np.mean(backgnd_seq['avoidance'].values),\
                      np.std(backgnd_seq['avoidance'].values)
    
    

    count = 0
    opt_df = pd.DataFrame()
    for sequence in mrna_df['sequence'] :
        message='at sequence :'+ str(count)
        functions.progress(count,mrna_df.shape[0],message)
        
        optimization = features.Optimize(sequence,cai_mean, cai_std,gc_cont_mean,\
                 gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,20)
        sequence_count = mrna_df.loc[mrna_df['sequence']\
                                                       == sequence].index[0]
        
        optimized_sequence = [optimization.simulated_anneal() \
                              for _ in range(2)]
        opt_df[sequence_count]=optimized_sequence
        count+=1
        message='at sequence :'+ str(count)
        functions.progress(count,mrna_df.shape[0],message)
        

        
    filename = mypath + o +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    opt_df.to_csv(filename,sep=',', encoding='utf-8', index=False)
    
    

    
                   
        
        

if __name__ == '__main__':
    m,o= check_arg(sys.argv[1:])
    main()
        
