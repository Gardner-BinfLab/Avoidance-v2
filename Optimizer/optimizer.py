import sys
import os
import argparse
import time
import pandas as pd
import numpy as np
import pickle
from os import path
from libs import functions,data,features
from multiprocessing import Pool, cpu_count
#from functools import partial


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Training script')
    parser.add_argument('-m', '--mrna',
                        type=valid_file,
                        help='mrna sequence from hmm in csv or fasta.',
                        required='True')
    parser.add_argument('-r', '--randomforest',
                        help='random forest model',
                        required='True')  
    parser.add_argument('-o', '--output',
                        help='Output file name.',
                        default = 'optimized_sequences')

    results = parser.parse_args(args)
    return (results.mrna,
            results.randomforest,
            results.output)


def all_features(dataframe):
    dataframe['analyze'] =  dataframe['sequence'].apply(lambda x:features.Analyze(x))
    dataframe['accs'] = dataframe['analyze'].apply(lambda x:x.access_calc())
    dataframe['sec_str'] = dataframe['analyze'].apply(lambda x:x.sec_str())
    dataframe['cai'] = dataframe['analyze'].apply(lambda x:x.cai_heg())
    dataframe['gc_cont'] = dataframe['analyze'].apply(lambda x:x.gc_cont())
    dataframe['avd'] = dataframe['analyze'].apply(lambda x:x.avoidance_opt())
    
    return dataframe

def parallelize_df(df, func):
    partitions = cpu_count()
    df_split = np.array_split(df, partitions)
    pool = Pool(partitions)
    results = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return results


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
        
    #first we find good sequences from the given list via random forest
    #then those good sequences can be optimized further by simulated annealing
    
    #for mrna scoring and selecting
    print('calculating features for given sequence list..', flush=True)
    mrna_df = parallelize_df(mrna_df,all_features)
    print('done!',flush=True)
    
    ##for random forest
    #load pre-trained rf
    print('loading random forest model..',flush=True)
    rf_model = pickle.load(open(r, 'rb'))
    print('done!',flush=True)
    
    
    print('selecting good sequences..', flush=True)
    mrna_df['rf_input'] = [mrna_df[['accs','sec_str','cai','gc_cont','avd']].values[x]\
                           for x in range(mrna_df.shape[0])]
    #we keep threshold of 0.9 i.e.. anything above or equal to 0.9 is 1, rest are 0
    mrna_df['rf_results'] = mrna_df['rf_input'].apply(lambda x: 1 if \
                                                      rf_model.predict_proba([x])[0][1] >=0.9 else 0)
    
    filename = mypath + 'mrna_analysis' +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    mrna_df.to_csv(filename,sep=',', encoding='utf-8', index=False)
    
    #pick those sequences with 1 from Random forest
    choosen_seq = mrna_df.loc[(mrna_df['rf_results'] == 1)].reset_index(drop=True)
    if choosen_seq.shape[0] == 0:
        print('failed to find any good sequences at this stage..',flush=True)
        print('optimizing all of provided sequence instead..', flush=True)
        choosen_seq = mrna_df
    else:
        filename = mypath + 'selected_sequences' +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
        choosen_seq.to_csv(filename,sep=',', encoding='utf-8', index=False)
        print('done!',flush=True)
    
    
    ##Now we generate 1000 random synonymous variants as a background
    #check if we already have background sequences and features
    if path.exists('background_sequences.csv'):
        print('previous background found! using it..', flush=True)
        backgnd_seq=pd.read_csv('background_sequences.csv')
    else:
        print('generating 1000 random synonymous sequences..', flush=True)
        backgnd_seq = functions.syn_background(mrna_df['sequence'][0],1000)

        #calculate features for background
        print('calculating features for background sequences..', flush=True)
        backgnd_seq = parallelize_df(backgnd_seq,all_features)
        backgnd_seq.to_csv('background_sequences.csv',index=False)
    
    
    #for z scores
    cai_mean, cai_std = np.mean(backgnd_seq['cai'].values),\
                        np.std(backgnd_seq['cai'].values)
    gc_cont_mean, gc_cont_std = np.mean(backgnd_seq['gc_cont'].values),\
                                np.std(backgnd_seq['gc_cont'].values)
    ss_mean, ss_std = np.mean(backgnd_seq['sec_str'].values),\
                      np.std(backgnd_seq['sec_str'].values)
    avd_mean, avd_std = np.mean(backgnd_seq['avd'].values),\
                      np.std(backgnd_seq['avd'].values)
    accs_mean, accs_std = np.mean(backgnd_seq['accs'].values),\
                      np.std(backgnd_seq['accs'].values) 
    
    
    #optimization begins
    print('optimization started..this may take a while..', flush=True)
    count = 0
    new_sequences = []
    for sequence in choosen_seq['sequence'] :
        message='at sequence :'+ str(count)
        functions.progress(count,choosen_seq.shape[0],message)
        
        optimization = features.Optimize(sequence,cai_mean, cai_std,gc_cont_mean,\
                 gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,accs_mean,accs_std,2000)
        
        
        pools = Pool(10)
        pool_results = []
        for result in pools.starmap(optimization.simulated_anneal,\
                                    [() for _ in range(10)]):
            pool_results.append(result)
        pools.close()
        pools.join()
        
        
        new_sequences.append(pool_results)
        count+=1
        message='at sequence :'+ str(count)
        functions.progress(count,choosen_seq.shape[0],message)
        
    
    #wew! finally we reached at the end
    print('optimization completed! exporting sequences...',flush=True)    
    final_sequences = pd.DataFrame({'sequence':[seq for sublist in new_sequences\
                                                for seq in sublist]})
    
    filename = mypath + o +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    final_sequences.to_csv(filename,sep=',', encoding='utf-8', index=False)
    print('full process was completed successfully!', flush = True)
    
                                    









if __name__ == '__main__':
    m,r,o= check_arg(sys.argv[1:])
    main()
        

