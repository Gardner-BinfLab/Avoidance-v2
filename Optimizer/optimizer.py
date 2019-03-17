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

#SA=Simulated Annealing


def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.csv', '.fasta','.fa'):
        raise argparse.ArgumentTypeError('File must have a csv or fasta extension')
    return param


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='optimizer script')
    parser.add_argument('-m', '--mrna',
                        type=valid_file,
                        help='mrna sequence from hmm in csv or fasta.',
                        required='True')
    parser.add_argument('-r', '--randomforest',
                        help='random forest model',
                        required='True')
    parser.add_argument('-u', '--utr5',
                        type=str,
                        help="5' utr (71 nt). default = pET")
    parser.add_argument('-s','--simanneal',
                        help='simulated annealing. pass this arg to do SA.',
                        action="store_true")
    parser.add_argument('-c','--count',
                        help='num of top sequences to pick for SA. default=10',
                        type=int,
                        default=10)
    parser.add_argument('-g','--gen',
                        help='num of times of SA per sequence. each SA gives one sequence. default 1',
                        type=int,
                        default=1)
    parser.add_argument('-n','--niter',
                        help='num itr in SA. default 200',
                        type=int,
                        default=200)
    parser.add_argument('-o', '--output',
                        help='output file name.',
                        default = 'optimized_sequences')


    results = parser.parse_args(args)
    return (results.mrna,
            results.randomforest,
            results.utr5,
            results.simanneal,
            results.count,
            results.gen,
            results.niter,
            results.output)


def all_features(dataframe):
    global utr_
    dataframe['analyze'] =  dataframe['sequence'].apply(lambda x:features.Analyze(x,utr=utr_))
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
    
    
    print('selecting',c,'good sequences..', flush=True)
    mrna_df['rf_input'] = [mrna_df[['accs','sec_str','cai','gc_cont','avd']].values[x]\
                           for x in range(mrna_df.shape[0])]
    #we keep threshold of 0.9 i.e.. anything above or equal to 0.9 is 1, rest are 0
    mrna_df['rf_scores'] = mrna_df['rf_input'].apply(lambda x:rf_model.predict_proba([x])[0][1])
    mrna_df_sorted = mrna_df.sort_values('rf_scores',ascending=False).reset_index(drop=True)  
    filename = mypath + o+ 'mrna_analysis' +'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv'
    mrna_df_sorted.to_csv(filename,sep=',', encoding='utf-8', index=False)
    
    #pick top 10 sequences from the analyzed sequences
    choosen_seq = mrna_df_sorted[:c]
    print('the max predicted probability of sequences being highly expressed is ',\
          choosen_seq['rf_scores'][0],flush=True)
    
    
    
    
    
    if s is True:
        
        ##Now we generate 1000 random synonymous variants as a background
        #check if we already have background sequences and features
        background_name = 'background_sequences_'+o+'.csv'
        if path.exists(background_name):
            print('previous background found! using it..', flush=True)
            backgnd_seq=pd.read_csv(background_name)
        else:
            print('generating 1000 random synonymous sequences as background..', flush=True)
            backgnd_seq = functions.syn_background(mrna_df['sequence'][0],1000)

            #calculate features for background
            print('calculating features for background sequences..', flush=True)
            backgnd_seq = parallelize_df(backgnd_seq,all_features)

            backgnd_seq.to_csv(background_name,index=False)


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
        for sequence in choosen_seq['sequence']:
            message='at sequence :'+ str(count)
            functions.progress(count,choosen_seq.shape[0],message)

            optimization = features.Optimize(sequence,cai_mean, cai_std,gc_cont_mean,\
                     gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,accs_mean,accs_std,500)

            pools = Pool(g)
            pool_results = []
            for result in pools.starmap(optimization.simulated_anneal,\
                                        [() for _ in range(g)]):
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
        
    else:
        pass
        
    print('full process was completed successfully!', flush = True)
    print(time.strftime("%Y%m%d-%H%M%S"))
                                    









if __name__ == '__main__':
    m,r,u,s,c,g,n,o= check_arg(sys.argv[1:])
    print('================================================')
    print(time.strftime("%Y%m%d-%H%M%S"))
    if u is None or len(u)!=71 :
        u = 'aggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacc'
    utr_ = u.lower()
    print("using ",utr_," as the 5' utr..", flush = True)
    if s is False:
        print('not doing simulated annealing..', flush = True)
    main()
        

