#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 20:17:27 2019

@author: bikash
"""

from subprocess import run, PIPE, DEVNULL
import tempfile
import os
import pandas as pd

class AnalyseTerminators:
    '''Analyses terminators for a dataframe of sequences using cmsearch
    '''


    def __init__(self, seq_df, cm):
        self.seq_df = seq_df
        self.cm = cm
        self.tempdir = os.path.join(tempfile.gettempdir(), 'cmsearch')
        self.tempfname = pd.util.testing.rands_array(15, 1)[0]
        self.cm_output = None
        self.results = None


    def make_rand_accs(self):
        '''Generate some random accession for sequences.
        '''
        self.seq_df['Accession'] = pd.util.testing.rands_array(10, len(self.seq_df))


    def make_tmp_dir(self):
        '''Make temporary directory for files from cmsearch.
        '''
        try:
            os.makedirs(self.tempdir)
        except FileExistsError:
            pass


    def dataframe_to_fasta(self):
        '''Export sequences to fasta. (Required for cmsearch input.)
        '''
        self.make_rand_accs()
        self.make_tmp_dir()
        file_contents = ''
        for ind, val in enumerate(self.seq_df.Accession):
            file_contents += '>'+val+'\n'+self.seq_df.Sequence[ind]+'\n'
        with open(self.tempdir+'/'+ self.tempfname + '.fa', 'w') as file_out:
            file_out.write(file_contents)


    def run_cmsearch(self):
        '''Run cmsearch
        '''
        self.dataframe_to_fasta()
        inp_f = self.tempdir + '/' + self.tempfname +'.fa'
        proc = run(['cmsearch', '--max', self.cm, inp_f],\
                   stdout=PIPE, stderr=DEVNULL, encoding='utf-8')
        os.remove(inp_f)
        self.cm_output = str(proc.stdout)


    def term_check(self):
        '''Parse results from cmsearch.
        It chops the table from cmsearch output and returns the dataframe
        with number of hits and E values.
        Note: No hits will have zero hits and high E value.
        '''
        if self.cm_output is None:
            self.run_cmsearch()

        tmp_res = self.cm_output.split('Hit scores:')[1].\
                    split('Hit alignments:')[0].split('\n')
        cmsearch_table = list(filter(None, tmp_res))[3:]

        accs = []
        e_val = []
        for _, val in enumerate(cmsearch_table):
            lst = list(filter(None, val.split(' ')))
            accs.append(lst[5])
            e_val.append(lst[2])


        cm_df_tmp = pd.DataFrame({'Accession':accs, 'E_val':e_val})

        #check for duplicates and group them
        cm_df = cm_df_tmp['E_val'].groupby(cm_df_tmp['Accession']).\
                apply(list).reset_index()

        #count number of hits.
        cm_df['Hits'] = cm_df['E_val'].apply(len)

        #find min eval (useful when sequence has multiple hits)
        cm_df['Min_E_val'] = cm_df['E_val'].apply(min)

        #merge with original dataframe
        final_results = pd.merge(self.seq_df, cm_df, on='Accession',\
                                 how='outer')

        #no hits are replaced by 0
        final_results['Hits'].fillna(0, inplace=True)

        #no hits are replaced by very high E-value
        final_results['Min_E_val'].fillna(10000, inplace=True)

        self.results = final_results.sort_values(['Hits', 'Min_E_val'],\
                                         ascending=[True, False]).\
                                         reset_index(drop=True)
        return self.results

