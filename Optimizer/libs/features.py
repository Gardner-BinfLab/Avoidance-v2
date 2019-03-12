
import os
import RNA
#from vienna import RNA
import tempfile
import numpy as np
import pandas as pd
from libs import data,functions
import subprocess
from subprocess import run,PIPE


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
        '''
        this does NOT cat ncrnas, rather spawn a process for each ncrna
        '''
        global ncrna
        sequence = self.sequence
        temp_df = pd.DataFrame({'input_mrna':['>input_mrna'+'\n' + sequence[:30]]})
        mrna_input = '\n'.join(temp_df['input_mrna'])
        try:
            ncrna.reset_index(level=0, inplace=True)
        except ValueError: #did we reset index already?
            pass
        ncrna['merge'] = ncrna['accession']+  ':break' +'\n' + ncrna['sequence'] 
        ncrna['input'] = ncrna['merge'] +'\n'+ mrna_input
        rnaup_res = functions.multiprocess_wrapper(functions.interaction_calc,\
                                                   ncrna['input'])
        avoidance = np.max(functions.rnaup_result_parser(rnaup_res)[0].values)
        return avoidance
    
    def avoidance_opt(self):
        '''
        this cats ncrnas and computes interaction for single mrna
        so is a single process
        '''
        global ncrna
        sequence = self.sequence
        #merge ncrnas to a single string
        try:
            ncrna.reset_index(level=0, inplace=True)
        except ValueError:
            pass
        ncrna['merge'] = ncrna['accession']+'\n' + ncrna['sequence'] +'\n'
        ncrna_input = ncrna['merge'].values.sum()
        mrna_input = '>input_mrna'+ ':break'+'\n' + sequence[:30] +'\n' +ncrna_input        
        rnaup_res = functions.interaction_calc(mrna_input)
    
        avoidance = np.max(functions.rnaup_result_parser(rnaup_res)[0].values)
        return avoidance
        
        
    def access_calc(self, length=30,
                    utr='ggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacat'\
                   ):


        tmp = os.path.join(tempfile.gettempdir(), 'plfold')
        try:
            os.makedirs(tmp)
        except FileExistsError:
            pass

        seq= utr.lower()+self.sequence.lower()
        proc = run(['RNAplfold', '-W 210', '-u 210', '-O'], \
                   stdout=PIPE,stderr=subprocess.DEVNULL,input=seq,\
                   cwd = tmp,encoding='utf-8') 

        open_en43 = pd.read_csv(tmp+'/plfold_openen',sep='\t',\
                               skiprows=2, header=None)[43][len(utr):len(utr)+length].sum()
        os.remove(tmp+'/plfold_openen')
        os.remove(tmp+'/plfold_dp.ps')
        return open_en43



        
        
    
class Optimize:
    '''does optimizations to a sequence
    '''
    def __init__(self,sequence,cai_mean, cai_std,gc_cont_mean,\
                 gc_cont_std,ss_mean, ss_std, avd_mean, avd_std,\
                 accs_mean,accs_std,niter=1000):
        
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
        self.accs_mean = accs_mean
        self.accs_std = accs_std
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

        avd_ = results.avoidance_opt()
        z_avd = Optimize.std_score(avd_, self.avd_mean, self.avd_std)
        
        accs_ = results.access_calc()
        z_accs = Optimize.std_score(accs_, self.accs_mean, self.accs_std)

        #total_z_score = z_cai - z_gc + z_ss + z_avd
        total_z_score = 2*z_cai - 0.1* z_gc + 0.5*z_ss + z_avd - z_accs

        return total_z_score
    
    
    def simulated_anneal(self):
        '''
        preforms a simulated annealing
        '''
        seq = self.sequence
        niter = self.niter
        temp = np.logspace(1,0.00001,niter)
        length = functions.sequence_length(seq)
        num_of_subst = [int((length-5)*np.exp(-_/int(niter/10))+5) \
                         for _ in range(niter)]
        scurr = seq
        sbest = seq
        for i in range(niter):
            T = temp[i]
            snew = functions.substitute_codon(sbest,num_of_subst[i])
            if self.cost_function(snew) >= self.cost_function(scurr):
                    scurr = snew
                    if self.cost_function(scurr)>=self.cost_function(sbest):
                        sbest = snew
            elif np.exp(-(self.cost_function(scurr)-self.cost_function(snew))/T)\
                            <= np.random.rand(1)[0]:
                scurr = snew
            text = '\ritr ' + '{:04d}'.format(i)
            message= f"\033[{'34m'}{text}\033[00m"
            print(message,end='\r')
        annealed_seq = sbest 
        return annealed_seq    



