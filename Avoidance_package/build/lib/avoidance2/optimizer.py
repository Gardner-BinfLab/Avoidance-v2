#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 13:47:59 2019

@author: bikash
"""

from functools import lru_cache
import numpy as np
from avoidance2 import functions, features




class Optimizer:
    '''Optimizes the sequence using simulated annealing
    '''


    def __init__(self, seq=None, niter=1000):
        self.seq = seq.upper()
        self.niter = niter
        self.back_feat = None
        self.back_mean_sd = None
        self.annealed_seq = None


    def background_features(self):
        '''Calculates features for background sequences.
        Currently, the background is random synonymous sequences.
        '''
        back_df = functions.syn_background(seq=self.seq, n=10)
        back_df_feat = features.AnalyzeDataFrameFeatures(input_file=back_df)
        back_mean_sd = back_df_feat.mean_and_sd()
        self.back_feat = back_df_feat.all_features
        self.back_mean_sd = back_mean_sd
        return back_mean_sd


    @lru_cache(maxsize=128, typed=True)
    def cost_function(self, new_seq=None):
        '''Calculates z scores of a given sequence based upon a group of
        sequences given to Optimizer class.
        '''
        if new_seq is None:
            seq = self.seq
        else:
            seq = new_seq
        ftrs = features.AnalyzeSequenceFeatures(seq=seq)
        if self.back_mean_sd is None:
            self.back_mean_sd = self.background_features()
        df = self.back_mean_sd
        df['weights'] = [1]*df.shape[0]
        df['seq'] = ''
        df.loc['cai', 'seq'] = ftrs.cai()
        df.loc['gc_cont', 'seq'] = ftrs.gc_cont()
        df.loc['sec_str', 'seq'] = ftrs.sec_str()
        df.loc['avd', 'seq'] = ftrs.avoidance_opt()
        df.loc['accs', 'seq'] = ftrs.access_calc()
        df['z_'] = (df['seq'] - df['mean'])/df['sd']
        df['weighted_z'] = df['z_'] * df['weights']
        total_z = df.loc[['cai', 'sec_str', 'avd', 'accs'], 'weighted_z'\
                         ].sum() - df.loc['gc_cont', 'weighted_z']
        return total_z


    def simulated_anneal(self):
        '''
        preforms a simulated annealing
        '''
        seq = self.seq
        niter = self.niter
        temp = np.geomspace(1, 0.00001, niter)
        length = functions.sequence_length(seq)
        num_of_subst = [int((length-5)*np.exp(-_/int(niter/10))+5) \
                         for _ in range(niter)]
        scurr = seq
        sbest = seq
        print("\nSimulated annealing progress:")
        print("iter\tcurrent cost\tbest cost")
        print("====\t============\t=========")
        for i in range(niter):
            t = temp[i]
            snew = functions.substitute_codon(sbest, num_of_subst[i])
            if self.cost_function(snew) >= self.cost_function(scurr):
                scurr = snew
                if self.cost_function(scurr) >= self.cost_function(sbest):
                    sbest = snew
            elif np.exp(-(self.cost_function(scurr)-self.cost_function(snew))\
                        /t) <= np.random.rand(1)[0]:
                scurr = snew
            print('\r%s\t %0.8f\t %0.8f\r'%(i, self.cost_function(scurr), \
                                             self.cost_function(sbest)), \
                                                end='\r')
        annealed_seq = sbest
        self.annealed_seq = annealed_seq
        return annealed_seq


