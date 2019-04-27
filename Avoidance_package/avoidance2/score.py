#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:06:17 2019

@author: bikash
"""

import os
import pickle
import warnings
import numpy as np
from avoidance2 import functions

class Score:
    '''
    Scores sequences based on their features using a random forest model.
    '''
    
    def __init__(self, input_file=None):
        self.input_file = functions.file_check(input_file)


    @staticmethod
    def validate_model(model):
        '''Checks if the given random forest classifier is usable.
        '''
        try:
            with open(model, 'rb') as file:
                rf = pickle.load(file)
        except TypeError:
            try:
                model.feature_importances_
                rf = model
            except AttributeError:
                raise NotImplementedError("Unknown model.")
        return rf
        
    
    def score_seq(self, rf=None):
        '''Scores sequence based on random forest classifier
        '''
        try:
            rf = self.validate_model(rf)
        except NotImplementedError:
            warnings.warn("Using default random forest classifier.")
            rf = self.validate_model(os.path.dirname(__file__)+\
                                              '/default_rf.pkl')
#        analysis_ = features.AnalyzeDataFrameFeatures(input_file=\
#                        self.input_file)
#        seq_features = analysis_.features()
        seq_features = self.input_file
        try:
            seq_features['rf_input'] = seq_features[['cai', 'gc_cont',\
                          'sec_str', 'avd', 'accs']].values.tolist()
        except KeyError:
            raise ValueError("Please calculate the sequences features first.")
        seq_features['rf_input'] = seq_features['rf_input'].apply(np.array)
        seq_features['rf_pred'] = seq_features['rf_input'].apply(\
                                lambda x: rf.predict_proba([x])[0][1])
        self.scored_seq = seq_features
        return seq_features
        

    def select_seq(self, num=None):
        '''Selects high scoring sequences based on classifier scores
        num specifies number of high scoring sequences to select
        '''
        try:
            num = int(num)
        except (ValueError, NameError, TypeError):
            num = 10
        seq_features = self.scored_seq
        seq_features.sort_values('rf_pred', inplace=True, ascending=False)
        seq_features = seq_features.reset_index(drop=True)
        selected_seq = seq_features[:num]
        print("The max scores of current selected seqences is {}.".\
              format(selected_seq['rf_pred'][0]))
        
        self.selected_seq = selected_seq
        return selected_seq


  