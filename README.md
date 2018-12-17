# Avoidance2.0
Alpha release 



# Install dependencies 
Requires Python3.

This command installs the required dependencies.

`pip3 install --user scikit-learn hmmlearn scipy numpy pandas progressbar2`


# Manual
`-h` gives possible options.

# Train
Train takes a list of sequences in CSV and builds the model. The built model is 
exported as a python pickle (.pkl) file in 'results/hmm' subdirectory.

Example use case with hidden states(`-n`)=2 and iterations (`-i`)=5

```sh
$ python3 train.py -f /home/bikash/Downloads/train_set_.csv -n 2 -i 5 -v
Reading input sequence file..
Seperating codons..
[=========================================================] 100% Time:  0:00:00

Initializing a model.
[=========================================================] 100% Time:  0:00:00
Training begins. This may take a while..
         1    -1802468.8592             +nan
         2    -1649991.9466     +152476.9126
         3    -1649984.8357          +7.1109
         4    -1649979.1202          +5.7155
         5    -1649973.6993          +5.4209
Training completed.
Model dumped.
```

# Scores
Scores outputs scores for given list of sequences given a model (`-m`) in .pkl
format. Note that scores != expression levels. They are just a measure of how 
likely the model is able to generate the given sequence. This is because the
model gives higher scores to soething which it knows and lower scores to
something which it doesn't. 

It's a good idea to find the baseline scores of training set to get an idea
of how much score is given by the model for 100% known sequences. If you check
scores of some sequences not used in training set but similar, then you'll get
scores closer to the baseline scores. If you give something different, the 
scores will be lower.

This means, the model can generate similar sequences in the neighbourhood of the
highly expressed sequences. However, similarity doesn't gurantee equal
expression and thus the generated sequences need to be checked. 

Also, note that for a model trained on highly expressed sequences,
upon increasing number of hidden states, it seems to start giving higher scores 
to "unknown" sequences, with higher expressions as well, while giving lower 
scores to those "unknown" sequences with poor  expressions.


Example use case:
```sh
$ python3 score.py -f /home/bikash/Downloads/train_set.pkl  -m example_model.pkl -v

Seperating codons..
[=========================================================] 100% Time:  0:00:00

Initializing a model.
[=========================================================] 100% Time:  0:00:00

Calculating Scores..
[=========================================================] 100% Time:  0:00:06
Exporting scores to csv..
```

# Generate
This generates sequences from given sequence. You have to specify the sequence,
model and number of sequence to generate. -v option recommended.


Example use case:
```sh
$ python3 seqgen.py -s 
agtgtttgtgtctgcaatcccaagtttgtttgcgctgaaatatgcgatgctcaatgttatgatctgcgtactaagccgcagatcatagtgggaact 
-m example_model.pkl -n 200 -v

Calculating Scores..
[=========================================================] 100% Time:  0:00:00
Exporting scores to csv..
```