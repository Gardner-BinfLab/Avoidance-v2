
import os
#import RNA
import numpy as np
import pandas as pd
from libs import data,functions
from multiprocessing import Pool
from subprocess import run, PIPE 

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
        global ncrna
        sequence = self.sequence
        mrna_input = '>input_sequence'+'\n'+sequence[:30]
        ncrna['input'] = ncrna.index+ ':break' + '\n'+ ncrna['sequence']\
                         + '\n'+ mrna_input
        ncrna['input_encoded'] = ncrna['input'].apply(lambda x: str.encode(x))
        rnaup_res = functions.multiprocess_wrapper(functions.interaction_calc,\
                                                   ncrna['input_encoded'])
        avoidance = functions.rnaup_result_parser(rnaup_res)[0]
        return avoidance
    
class Optimize:
    '''does optimizations to a sequence
    '''
    def __init__(self,sequence,background_df,cai_mean, cai_std,gc_cont_mean,\
                 gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,niter=1000)
    self.sequence = sequence.lower()
    self.background_df = background_df
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
    
    
    def cost_function(self):
        
        @staticmethod
        def std_score(x,mu,sigma):
            z = (x-mu)/sigma
            return z
        
        sequence = self.sequence
        results = Analyze(sequence)
        
        cai_ = results.cai()
        z_cai = std_score(cai_, cai_mean, cai_std)
        
        gc_ = results.gc_cont()
        z_gc = std_score(gc_, gc_cont_mean, gc_cont_std)


        ss_ = results.sec_str()
        z_ss = std_score(ss_, ss_mean, ss_std)

        avd_ = results.avoidance()
        z_avd = std_score(avd_, avd_mean, avd_std)

        total_z_score = z_cai + z_gc + z_ss + z_avd

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
            if cost_function(snew) >= cost_function(scurr):
                    scurr = snew
                    if cost_function(scurr)>=cost_function(sbest):
                        sbest = snew
            elif np.exp(-(cost_function(scurr)-cost_function(snew))/T)\
                            <= np.random.rand(1)[0]:
                scurr = snew
            functions.progress(i+1,niter)
        annealed_seq = sbest 
        return annealed_seq    



