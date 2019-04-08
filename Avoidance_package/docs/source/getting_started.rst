***************
Getting started
***************

Installation
============

**Prerequisite**

You need to install Python 3.6 and ViennaRNA. Installation guide for ViennaRNA
is at `their website. <https://www.tbi.univie.ac.at/RNA/>`_ Please make sure 
that you install 'Python 3 bindings' as well. 
Once those steps are completed, you can install our software with this simple 
command.

.. code::

    pip3 install --user avoidance2

Ths will automatically install the dependencies `Pandas 0.23.4` and 
`scikit-learn 0.20.2`. 

Sampling your first sequences
================================

To sample synonymous sequences, first, you need to import MarkovModel from the 
package and initialize it.

.. code:: 

    from avoidance2.markov import MarkovModel
    my_model = MarkovModel() #initialize

Now, given an input sequence, we can sample any number of synonymous sequences.
In this example, we sample 10 synonymous sequences.

.. code:: 

    my_model.generate(sequence='ATGCAGCGAGCGAGCGAGCGAGCAGC', num_seq=10)
    |██████████████████████████████████████| 100% (10/10) Synonymous sequences 
    Completed!
    
                       sequence
    0  ATGCAACGTGCATCTGAACGTAGC
    1  ATGCAACGCGCTAGTGAACGCAGC
    2  ATGCAACGTGCAAGTGAACGTAGC
    3  ATGCAACGTGCTTCTGAGCGTAGC
    4  ATGCAGCGCGCTTCTGAACGCAGC
    5  ATGCAACGAGCGAGTGAGCGAAGC
    6  ATGCAACGTGCTTCGGAACGTAGC
    7  ATGCAAAGAGCTTCTGAAAGAAGC
    8  ATGCAACGTGCGTCTGAGCGTAGC
    9  ATGCAACGAGCGTCTGAGCGAAGC
    
This returns a dataframe of synonymous sequences, which can be easily exported
to a csv file.

.. code:: 

    import pandas as pd
    test_seq = 'ATGCAGCGAGCGAGCGAGCGAGCAGC
    syn_seq = my_model.generate(sequence=test_seq, num_seq=10)
    syn_seq.to_csv('syn_seq_test.csv', index=None)

This example uses the default model provided with the package to generate
sequences. You can build your own model by using your own sequences. In this
example, we train the model for upto 50 codons (0 to 49).

.. code:: 

    from avoidance2.markov import MarkovModel
    new_mm = MarkovModel(input_file='markov_test.csv', length=50)
    mm_built = new_mm.train()
    |██████████████████████████████████████| 100% (49/49) Training.. 
    Completed!


This model can be used to generate sequences.

.. code:: 

    new_mm.generate('ATGAGCGAGCGAGGCGAGCGACGAGCGAGCGACGG', num_seq=10)
    :
    UserWarning: Out of 10 generated sequences, 2 duplicate sequences were 
    removed.
                                sequence
    0  ATGTCAGAGCGCGGCGAGCGCCGCGCTTCACGG
    1  ATGTCCGAGCGCGGCGAGCGCCGCGCCTCCCGG
    2  ATGTCCGAGCGGGGTGAGCGGCGGGCGTCCCGG
    3  ATGTCCGAGAGGGGCGAGAGGAGGGCGTCCCGG
    4  ATGTCCGAGAGAGGCGAGAGAAGAGCGTCCCGG
    5  ATGTCCGAGCGGGGCGAGCGGCGGGCATCCCGG
    6  ATGTCCGAACGAGGCGAACGACGAGCGTCCCGG
    7  ATGTCCGAGAGGGGGGAGAGGAGGGCGTCCCGG


    
The built model can be exported for future use using 'pickle'.

.. code:: 

    import pickle
    with open('markov.pkl', 'wb') as model:
        pickle.dump(mm_built, model)


Reusing this exported model can be done as follows:

.. code:: 

    new_mm.generate('ATGAGCGAGCGAGGCGAGCGACGAGCGAGCGACGG', model='markov.pkl',
     num_seq=10)
    
  

Computing features for your sequence
=====================================

We currently compute codon adaptation index (cai), G+C content (gc_cont), 
secondary structure (sec_str), avoidance (avd) and accessibility (accs).
To compute features for your sequence, you can do the following:

.. code:: 

    from avoidance2.features import AnalyzeSequenceFeatures
    sequence = 'ATGAGAGGCGAGCAGAGCGAGCGACGACGACGACGCAGAGCG'
    analysis_ = AnalyzeSequenceFeatures(sequence)

Your desired features can be found by calling the analysis object that we just 
created.

.. code:: 

    analysis_.cai() #Codon adaptation index
    : 0.04187077422415583
    analysis_.gc_cont() #G+C content
    : 0.5714285714285714
    analysis_.sec_str() #Secondary structure
    : -7.199999809265137
    analysis_.avoidance() #Avoidance
    : -4.82
    analysis_.access_calc() #Accessibility
    : 287.058339
    
    
Building a random forest model
==============================

Once the features of your sequences are calculated, we can use it to buld a 
random forest. The random forest classifier can be used to rank and predict 
the probability of being expressed.

.. code:: 

    from avoidance2.classifier import RandomForest
    rfc = RandomForest(input_file='rf_test.csv') #initialize the model
    my_rand_f = rfc.build_model() #build it
    :
    Random forest summary:
    n_est    accuracy        AUC     
    =====    ========        ====    
    50       0.650           0.702  
    
    Gini features importance:
    ==========================
    ['avd', 'cai', 'sec_str', 'gc_cont', 'accs'] 
     [0.1609 0.1769 0.1955 0.2324 0.2342]


The built model can be exported for future use using 'pickle'.

.. code:: 

    import pickle
    with open('new_rfc.pkl', 'wb') as model:
        pickle.dump(my_rand_f, model)
    



Scoring the sequence using random forest
===========================================

Sequences whose features are already calculated can be scored by using a 
random forest. This example uses the 'my_rand_f' model built above.

.. code:: 

    from avoidance2.score import Score
    sc = Score(input_file='rf_test.csv') #initalize score
    sc.score_seq(rf=my_rand_f) #scoring using the random forest
    :
            accession   ...   rf_pred
    0   >PaCD00423476   ...      0.26
    1   >PaCD00651452   ...      0.92
    2   >EcCD00344566   ...      0.96
    3   >HmCD00598756   ...      0.32
    4   >AtCD00589313   ...      0.78
    5   >HmCD00536337   ...      0.96
    ... ...             ...      ...
    

Once the sequences are scored, you can just pick the top 'n' number of
sequences, if you wish.

.. code:: 

    sc.select_seq(5) #picking 5 top sequences
    :The max scores of current selected seqences is 1.0.
           accession   ...   rf_pred
    0  >BsCD00605596   ...       1.0
    1  >VpCD00334889   ...       1.0
    2  >SoCD00338340   ...       1.0
    3  >BtCD00339745   ...       1.0
    4  >MmCD00540649   ...       1.0
    
This is a dataframe which can be exported easily as above.


Optimizing sequences
====================

If your sequences are scoring lower, you can optimize them further. The 
optimization by simulated annealing. We will take the sequence '>PaCD00423476'
which scored 0.26 in our model.

.. code:: 

    from avoidance2.optimization import Optimizer
    poor_seq = 'ATGGCCCATGCCGCGGCGGTCGAGGAAAACCAGTTGCAGTCGCTGAAGGACGT...'
    op = Optimizer(poor_seq, my_rand_f) #using the previous random forest model
    optimized_seq = op.simulated_anneal(itr=10) #10 iterations

    Simulated annealing progress:
    iter    current cost    best cost
    ====    ============    =========
    9        1.00000000      1.00000000
    
Our optimized sequence has scored 1 on the random forest!

