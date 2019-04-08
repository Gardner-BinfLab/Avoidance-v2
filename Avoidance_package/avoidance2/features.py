#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 18:51:16 2019

@author: bikash
"""

import os
import warnings
from functools import partial
import tempfile
from subprocess import run, PIPE, DEVNULL
import numpy as np
import pandas as pd
import RNA
from avoidance2 import data, functions, config



class AnalyzeSequenceFeatures():
    '''Analyze sequence features
    '''

    def __init__(self, seq):
        self.seq = seq.upper()



    def cai(self):
        '''Codon adaptation index
        '''
        given_seq = functions.splitter(self.seq[3:-3])
        excluded_codons = ['ATG', 'TGG']
        codon_list = [codon for codon in given_seq if codon not in \
                      excluded_codons]
        if any(codon in codon_list for codon in data.STOP_CODONS):
            raise KeyError("Given sequence has stop codons inside coding"
                           " region.")
        try:
            cai_values = [np.log(data.cai_table[codon]) for codon\
                          in codon_list]
            score = np.exp(np.mean(cai_values))
        except Exception as exc:
            raise KeyError("Either the given sequence is malformed or"
                           "the lookup table is corrupted.") from exc

        return score


    def gc_cont(self):
        '''G+C content
        '''
        g_count = self.seq[3:-3].count('G')
        c_count = self.seq[3:-3].count('C')
        gc_cont = (g_count + c_count)/len(self.seq)
        return gc_cont


    def sec_str(self, utr=None, positions=None):
        '''Secondary structure of specified position in sequence.
        Default position is -30 to 30 with 0 as start codon
        '''
        if utr is None:
            utr = 'GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
            if config.UTR_WARN:
                warnings.warn("Using default T7 promoter sequence from pET "
                              "vector as UTR.")
                config.UTR_WARN = False
        else:
            utr = utr.upper()
        if positions is None:
            if config.POSITION_WARN:
                warnings.warn("Looking for secondary structures in 30 codons "
                              "before and after start codon.")
                config.POSITION_WARN = False
            positions = (-30, 30)
        sequence = self.seq
        utr = utr[positions[0]:]
        seq = sequence[:positions[1]]
        total_sequence = utr+seq
        mfe = RNA.fold(total_sequence)[1]
        return mfe


    def avoidance(self, ncrna_file=None):
        '''
        This cats ncrnas and computes interaction for single mrna
        so is a single process
        '''
        if ncrna_file is None:
            if config.NCRNA_WARN:
                warnings.warn("Using default 5S RNA.")
                config.NCRNA_WARN = False

            ncrna_file = os.path.dirname(__file__) + '/ncrna.fa'
        ncrna = functions.file_check(ncrna_file)
        sequence = self.seq
        #merge ncrnas to a single string
        try:
            ncrna.reset_index(level=0, inplace=True)
        except ValueError:
            pass
        ncrna['merge'] = ncrna['accession']+'\n' + \
                                ncrna['sequence'] +'\n'
        ncrna_input = ncrna['merge'].values.sum()
        mrna_input = '>input_mrna'+ ':break'+'\n' + sequence[:30] +'\n' + \
                        ncrna_input
        rnaup_res = functions.interaction_calc(mrna_input)
        avoidance = np.min(functions.rnaup_result_parser(rnaup_res)[0].values)
        return avoidance


    def access_calc(self, length=None, utr=None):
        '''Sequence accessibility
        '''
        tmp = os.path.join(tempfile.gettempdir(), 'plfold')
        try:
            os.makedirs(tmp)
        except FileExistsError:
            pass
        if length is None:
            length = 30
        if utr is None:
            utr = 'GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
            if config.UTR_WARN:
                warnings.warn("Using default T7 promoter sequence from pET "
                              "vector as UTR.")
                config.UTR_WARN = False
        else:
            utr = utr.upper()
        sequence = utr + self.seq
        seq_accession, rand_string = functions.accession_gen()
        input_seq = seq_accession + sequence
        run(['RNAplfold', '-W 210', '-u 210', '-O'], \
                   stdout=PIPE, stderr=DEVNULL, input=input_seq, cwd=tmp, \
                    encoding='utf-8')
        out1 = '/' + rand_string + '_openen'
        out2 = '/' + rand_string + '_dp.ps'
        open_en43 = pd.read_csv(tmp+out1, sep='\t', skiprows=2, header=None)\
                    [43][len(utr):len(utr)+length].sum()

        os.remove(tmp+out1)
        os.remove(tmp+out2)
        return open_en43



class AnalyzeDataFrameFeatures():
    '''Calculates features for sequences in a dataframe
    '''


    def __init__(self, input_file=None, utr=None, positions=None, \
                 ncrna_file=None, length=None, multiprocess=True):
        self.input_file = functions.file_check(input_file)
        self.multiprocess = multiprocess
        self.utr = utr
        self.positions = positions
        self.ncrna_file = ncrna_file
        self.length = length
        self.all_features = None
        self.mean_sd = None


    def features(self):
        '''Calculates all features for given sequences
        '''
        input_df = self.input_file
        ftrs_func = partial(functions.sequence_df_features, self.utr, \
                            self.positions, self.ncrna_file, self.length)
        if self.multiprocess:
            results = functions.parallelize_df(input_df, ftrs_func)
        else:
            results = functions.sequence_df_features(self.utr, \
                    self.positions, self.ncrna_file, self.length, input_df)
        self.all_features = results
        return results



    def mean_and_sd(self):
        '''Calculates mean and standard deviation of features in a dataframe
        '''
        if self.all_features is None:
            self.all_features = self.features()
        df_with_features = self.all_features
        mean_sd = pd.DataFrame({'mean':df_with_features.mean(), \
                             'sd':df_with_features.std(ddof=0)})
        #ddof is Bessel's correction. We turn that off
        self.mean_sd = mean_sd
        return mean_sd

