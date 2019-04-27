#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 16:55:27 2019

@author: bikash
"""

from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from avoidance2 import functions, config

class RandomForest:
    '''Random forest classifier
    '''

    def __init__(self, input_file=None, n_est=None):
        self.input_file = functions.file_check(input_file)
        self.n_est = n_est
        if self.n_est is None:
            self.n_est = 50
        self.rfc = None


    @staticmethod
    def auc(X_test, y_test, model):
        '''Area under a ROC curve
        '''
        probs = model.predict_proba(X_test)
        preds = probs[:, 1]
        fpr, tpr = metrics.roc_curve(y_test, preds)[:2] #3rd ->thresholds
        auc = metrics.auc(fpr, tpr)
        return auc


    def build_model(self):
        '''Builds a random forest classifier.
        '''
        features = ['cai', 'gc_cont', 'sec_str', 'avd', 'accs']
        try:
            X = self.input_file[features].values
        except KeyError:
            raise ValueError("The given file doesn't have required features."
                             " Please calculate the features first.")
        try:
            y = self.input_file['label'].values
        except KeyError:
            raise ValueError("The given file has no binary classification "
                             "label. Please ensure your file has a column "
                             " 'label' with binary values (1 or 0)")
        X_train, X_test, y_train, y_test = train_test_split(X, y, \
                                            test_size=0.2, \
                                            random_state=config.RANDOM)
        clf = RandomForestClassifier(n_estimators=self.n_est, \
                                          random_state=config.RANDOM)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy = metrics.accuracy_score(y_test, y_pred)
        auc = self.auc(X_test, y_test, clf)
        imp_feat = [features[i] for i in np.argsort(clf.feature_importances_)]
        feat_gini_sorted = np.round(np.sort(clf.feature_importances_), 4)
        print("\nRandom forest summary:")
        print("n_est\t accuracy\t AUC\t ")
        print("=====\t ========\t ====\t ")
        print('\r%d\t %0.3f\t\t %0.3f\t\n'%(self.n_est, accuracy, auc), \
              end='\r')                                     
        print("\nGini features importance:")
        print("==========================")
        print(imp_feat, "\n", feat_gini_sorted)
        self.rfc = clf
        return clf
        