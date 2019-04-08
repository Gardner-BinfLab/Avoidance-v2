#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 12:49:12 2019

@author: bikash
"""

import os
import pickle
import warnings
import pandas as pd
import numpy as np
from avoidance2 import functions, data, config


class MarkovModel:
    '''This is used to train a Markov model using the input sequence.
    By default, we train for length 32 codons. The start codon and stop codon
    are stripped by default, but if your sequences are already free from them,
    you can use nostartstop=False.
    '''

    def __init__(self, input_file=None, nostartstop=True, length=32):
        try:
            self.input_file = functions.file_check(input_file)
        except NotImplementedError:
            self.input_file = None
        self.nostartstop = nostartstop
        self.length = length
        self.model = None
        self.synonymous_sequences = None


    def codons_to_df(self):
        '''
        Returns a dataframe of codons in given sequences. position is column
        and sequence number is index.
        '''
        codon_df = pd.DataFrame()
        training_seq_df = self.input_file
        if self.input_file is None:
            raise ValueError("Input file is not provided. Nothing to do!")
        if self.nostartstop is not True:
            warnings.warn("Start/stop codons not stripped. To reduce the "
                          "possibility of future errors, please ensure "
                          "that your sequences are already free from start"
                          " and stop codons or run the program again with "
                          "this flag set to True.")
        codon_df['sequence'] = training_seq_df['sequence'].values
        codon_df['codons'] = codon_df['sequence'].apply(functions.splitter)
        if self.nostartstop is True:
            codon_df['codons'] = codon_df['codons'].\
                                    apply(lambda x: x[1:-1])
        total_seq = codon_df.shape[0]
        codon_df['codons'] = codon_df['codons'].\
                                apply(lambda x: np.nan if any(codon in \
                                    x[3:-3] for codon in data.STOP_CODONS)\
                                    else x)
        codon_df = codon_df.dropna()
        remained_seq = codon_df.shape[0]
        if total_seq != remained_seq:
            warnings.warn("{} sequences were removed because we found stop"
                          "codons inside the coding region.".format(total_seq-\
                                                             remained_seq))
        splitted_df = pd.DataFrame(codon_df['codons'].values.tolist(), \
                                   index=codon_df.index)

        return splitted_df


    def train(self):
        '''
        Uses the codon dataframe to calculate probabilities
        codon name is the index and column number is the position of codon
        along the sequence.
        '''

        codon_df = self.codons_to_df()
        length = self.length
        if codon_df.shape[1] < self.length:
            length = codon_df.shape[1]
        try:
            codon_df = codon_df[list(range(int(length)))]
        except ValueError:
            pass
        codons = [key for key, value in data.codon2aa.items()]
        prob_df = pd.DataFrame(index=codons)
        for i in codon_df.columns:
            prob_df[i] = np.nan*len(prob_df)
            probs = codon_df.groupby(i).size().\
                                          div(codon_df[i].count())
            for item in codon_df[i]:
                if str(item).lower() not in ['nan', 'none']:
                    prob_df[i].loc[item] = probs.loc[item]
            functions.progress(i, len(codon_df.columns)-1, \
                               message="Training..")

        #fill nans with codon probability in entire dataset
        counts = codon_df.apply(pd.value_counts)
        total = counts.fillna(0).values.sum()
        prob_df['codon_prob'] = np.nan
        #fill with codon probs
        for i in range(len(counts.index)):
            prob_df['codon_prob'][counts.index[i]] = counts.\
                                                        loc[counts.index[i]].\
                                                        sum()/total
       #mask NaNs with codon probs
        for i in range(prob_df.shape[1]-1):
            prob_df[i].fillna(prob_df['codon_prob'], inplace=True)
        #normalize
        for j in range(prob_df.shape[1]-1):
            prob_df[j] = prob_df[j]/prob_df[j].sum()
        self.model = prob_df
        return prob_df


    @staticmethod
    def validate_model(model):
        '''Checks if the given Markov model is usable.
        '''
        codons = [k for k, v in data.codon2aa.items()]
        try:
            with open(model, 'rb') as file:
                prob_df = pickle.load(file)
            try:
                all(prob_df.index == codons)
            except Exception as exc:
                raise LookupError("Required codons were not found in the "
                                  "given model. Probably it is an older model"
                                  ". Please build and export a fresh model "
                                  "again.") from exc
        except TypeError:
            try:
                test_pos = np.random.randint(0,61)
                if all(model.iloc[test_pos].values == \
                       model.loc[codons[test_pos]].values[0]):
                    prob_df = model
            except (AttributeError, NotImplementedError):
                raise NotImplementedError("Unknown model.")
        return prob_df


    @staticmethod
    def get_probs(possible_codons, index, prob_df):
        '''Tries to extract probabilities of possible codons for substitution
        at given index(position along the sequence) using a probability table.
        First a lookup at that position is done, if that fails, we lookup to
        'codon_prob' column which is the overall probability of that codon in
        the training sequence (priors). If that fails, we return a random probability
        using 1D Dirichlet distribution.
        '''
        try:
            probs = prob_df.loc[possible_codons.substit[index], index]/\
                    prob_df.loc[possible_codons.substit[index], index].sum()
        except (KeyError, RuntimeWarning, IndexError):
            probs = prob_df.loc[possible_codons.substit[index], 'codon_prob']\
                                /prob_df.loc[possible_codons.substit[index], \
                                            'codon_prob'].sum()
        if probs.dropna().empty:
            probs = np.random.dirichlet\
                    (np.ones(len(possible_codons.substit[index])))
        return probs


    def generate(self, sequence, model=None,\
                         num_seq=10000):
        '''Generate synonymous sequences to a given sequence using a Markov
        model
        '''
        try:
            prob_df = self.validate_model(model)
        except NotImplementedError:
            if self.model is None:
                if config.DEFAULT_MARKOV:
                    warnings.warn("You haven't trained the model yet, nor "
                                  "supplied the trained model. So we will use "
                                  "the default model to generate sequences.")
                    config.DEFAULT_MARKOV = False
                prob_df = self.validate_model(os.path.dirname(__file__)+\
                                              '/default_markov.pkl')
            else:
                prob_df = self.model
        seq = sequence.upper()
        start = seq[:3]
        stop = seq[-3:]
        codons = functions.splitter(seq[3:-3])
        try:
            amino_acids = [data.codon2aa[c] for c in codons]
        except Exception as exc:
            raise ValueError("Stop codons were found inside the sequence.") from\
                exc

        possible_codons = pd.DataFrame({'substit':[data.aa2codon[aa] for aa\
                                                   in amino_acids]})
        subst_df = pd.DataFrame(columns=amino_acids)
        for index, a_a in enumerate(subst_df.columns):
            probs = self.get_probs(possible_codons, index, prob_df)
            try:
                subst_df[a_a] = np.random.choice(possible_codons.\
                        substit[index], num_seq, p=probs)
            except RuntimeWarning as exc:
                raise ValueError("We've encountered a serious problem. Please "
                                 "report this issue to the developers with "
                                 "the training sequences and sequence you are "
                                 "trying to generate synonymous sequences "
                                 "with.") from exc
            functions.progress(int((index+1)*num_seq/subst_df.shape[1]),\
                               num_seq, message="Synonymous sequences")


        subst_df['sequence'] = pd.Series(subst_df.values.tolist()).map(lambda\
                                 x: ''.join(map(str, x))).apply(lambda x: \
                                           start + x + stop)
        total_seq = subst_df.shape[0]
        final_df = pd.DataFrame({'sequence':subst_df['sequence']}).\
                    drop_duplicates().reset_index(drop=True)
        remained_seq = final_df.shape[0]
        if total_seq != remained_seq:
            warnings.warn("Out of {} generated sequences, {} duplicate "
                          "sequences were removed."\
                          .format(total_seq, total_seq-remained_seq))
        self.synonymous_sequences = final_df
        return final_df

