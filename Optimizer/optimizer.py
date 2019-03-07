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
                        default = 'optimized sequences')

    results = parser.parse_args(args)
    return (results.background,
            results.mrna,
            results.output)






def main():
    np.random.seed(12345)
    
    mypath = os.path.join(os.getcwd(),'results','model','')
    if os.path.exists(mypath)==True:
        pass
    else:
        os.makedirs(os.path.join(os.getcwd(),'results','model',''))

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
    backgnd_seq['features'] = backgnd_seq['sequence'].apply(lambda x:\
                                                            features.Analyze(x))
    backgnd_seq['cai'] = backgnd_seq['features'].apply(lambda x : x.cai())
    backgnd_seq['gc_cont'] = backgnd_seq['features'].apply(lambda x : x.gc_cont())
    backgnd_seq['sec_str'] = backgnd_seq['features'].apply(lambda x : x.sec_str())
    backgnd_seq['avoidance'] = backgnd_seq['features'].apply(lambda x :\
                                                             x.avoidance())
    
    #for z scores
    cai_mean, cai_std = np.mean(backgnd_seq['cai'].values),\
                        np.std(backgnd_seq['cai'].values)
    gc_cont_mean, gc_cont_std = np.mean(backgnd_seq['gc_cont'].values),\
                                np.std(backgnd_seq['gc_cont'].values)
    ss_mean, ss_std = np.mean(backgnd_seq['sec_str'].values),\
                      np.std(backgnd_seq['sec_str'].values)
    avd_mean, avd_std = np.mean(backgnd_seq['avoidance'].values),\
                      np.std(backgnd_seq['avoidance'].values)
    
    
    for sequences in mrna_df['sequence'].iterrows():
        optimization = features.Optimize(sequence)
        optimizied_sequence = optimization.simulated_anneal()
        ##to_do
    

        
        

if __name__ == '__main__':
    b,m,o= check_arg(sys.argv[1:])
    main()
        