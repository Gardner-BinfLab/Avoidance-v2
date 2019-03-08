
import os
import RNA
import numpy as np
import pandas as pd
from libs import data,functions


script_dir = os.path.dirname(__file__)
ncrna_path = os.path.join(script_dir, './ncrna.fa')
ncrna = functions.read_fasta(ncrna_path)


class Analyze():
    '''analyze sequence features
    '''
    def __init__(self,sequence,positions=(-30,30),\
                 utr='ggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacat'):
        self.sequence = sequence.lower()
        self.utr = utr
        self.positions = positions
        
        
        
    def cai(self):
        seq = self.sequence
        given_seq = functions.splitter(seq,len(seq))
        try:
            cai_values = [np.log(data.cai_table[codon]) for codon\
                          in given_seq]
            score = np.exp(np.mean(cai_values))
        except KeyError:
            print('strange sequence or corrupted cai table!')
            return 0
        return score

    
    def gc_cont(self):
        seq = self.sequence
        g_count = seq.count('g')
        c_count = seq.count('c')
        gc_cont = (g_count + c_count)/len(seq)
        return gc_cont

    
    def sec_str(self):
        utr = self.utr
        positions = self.positions
        sequence = self.sequence
        utr = utr[:positions[0]]
        seq = sequence[:positions[1]]
        total_sequence = utr+seq
        ss,mfe = RNA.fold(total_sequence)
        return mfe
    
    
    def avoidance(self):
        '''for single mrna vs many ncrnas!
        '''
        global ncrna
        sequence = self.sequence
        temp_df = pd.DataFrame({'input_mrna':['>input_mrna'+'\n' + sequence[:30]]})
        mrna_input = '\n'.join(temp_df['input_mrna'])
        ncrna['input'] = ncrna.index+ ':break' + '\n'+ ncrna['sequence'] + \
                        '\n'+ mrna_input
        ncrna['input_encoded'] = ncrna['input'].apply(lambda x: str.encode(x))
        rnaup_res = functions.multiprocess_wrapper(functions.interaction_calc,\
                                                   ncrna['input_encoded'])
        avoidance = np.max(functions.rnaup_result_parser(rnaup_res)[0].values)
        return avoidance
    
class Optimize:
    '''does optimizations to a sequence
    '''
    def __init__(self,sequence,cai_mean, cai_std,gc_cont_mean,\
                 gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,niter=1000):
        
        self.sequence = sequence.lower()
        self.cai_mean = cai_mean
        self.cai_std = cai_std
        self.gc_cont_mean = gc_cont_mean
        self.gc_cont_std = gc_cont_std
        self.ss_mean = ss_mean
        self.ss_std = ss_std
        self.avd_mean = avd_mean
        self.avd_mean = avd_mean
        self.avd_std = avd_std
        self.niter = niter
    
    
    @staticmethod
    def std_score(x,mu,sigma):
        z = (x-mu)/sigma
        return z
    
    
    
    def cost_function(self,new_seq=None):

        
        if new_seq is None:
            sequence = self.sequence
        else:
            sequence = new_seq
        results = Analyze(sequence)
        
        cai_ = results.cai()
        z_cai = Optimize.std_score(cai_, self.cai_mean, self.cai_std)
        
        gc_ = results.gc_cont()
        z_gc = Optimize.std_score(gc_, self.gc_cont_mean, self.gc_cont_std)


        ss_ = results.sec_str()
        z_ss = Optimize.std_score(ss_, self.ss_mean, self.ss_std)

        avd_ = results.avoidance()
        z_avd = Optimize.std_score(avd_, self.avd_mean, self.avd_std)

        total_z_score = z_cai - z_gc + z_ss + z_avd

        return total_z_score
    
    
    def simulated_anneal(self):
        '''
        preforms a simulated annealing
        '''
        seq = self.sequence
        niter = self.niter
        temp = np.linspace(1,0.001,niter)
        scurr = seq
        sbest = seq
        for i in range(niter):
            T = temp[i]
            snew = functions.substitute_codon(sbest)
            if self.cost_function(snew) >= self.cost_function(scurr):
                    scurr = snew
                    if self.cost_function(scurr)>=self.cost_function(sbest):
                        sbest = snew
            elif np.exp(-(self.cost_function(scurr)-self.cost_function(snew))/T)\
                            <= np.random.rand(1)[0]:
                scurr = snew
            functions.progress(i+1,niter)
        annealed_seq = sbest 
        return annealed_seq    



