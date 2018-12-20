# Avoidance2.0
Alpha release.
All changes will be made in testing branch and merged to master (this branch).



# Install dependencies 
Requires ~~`Python3.6+`~~ `Python3.4+`. (commits ad6588c, c90ba69, ad6588c)

This command installs the required dependencies.

`pip3 install --user scikit-learn hmmlearn scipy numpy pandas progressbar2`


# Manual
  - `-h` gives possible options.
  - `-v` highly recommended (for checking progress at some lengthy steps)
  - An example model is included and some training sequences are in test folder.


# Train
Train takes a list of sequences in CSV and builds the model. The built model is 
exported as a python pickle (.pkl) file in 'results/hmm' subdirectory.

Example use case with hidden states(`-n`)=2 and iterations (`-i`)=5

```console
$ python3 train.py -f /home/bikash/Downloads/train_set_.csv -n 2 -i 5 -v
Reading input sequence file..
Separating codons..
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
similar is the given sequence to the training sequence. This is because the
model gives higher scores to something which it knows and lower scores to
something which it doesn't. 

It's a good idea to find the baseline scores of training set to get an idea
of how much score is given by the model for 100% known sequences. If you check
scores of some sequences not used in training set but similar, then you'll get
scores closer to the baseline scores. If you give something different, the 
scores will be lower.

This means, the model can generate similar sequences in the neighbourhood of the
training sequences. However, similarity doesn't gurantee equal
expression and thus the generated sequences need to be checked. 


Example use case:
```console
$ python3 score.py -f /home/bikash/Downloads/train_set.pkl  -m example_model.pkl -v

Seperating codons..
[=========================================================] 100% Time:  0:00:00

Initializing a model.
[=========================================================] 100% Time:  0:00:00

Calculating Scores..
[=========================================================] 100% Time:  0:00:06
Exporting scores to csv..
```

Also, note that for a model trained on highly expressed sequences,
upon increasing number of hidden states, it seems to start giving higher scores 
to "unknown" sequences, with higher expressions as well, while giving lower 
scores to those "unknown" sequences with poor  expressions.

In the figure below red, teal and blue are "good","average" and "bad" sequences
from the training set. The model was trained on "good" ones and it does not _know_ 
the "average" and "bad" sequences. However, since they are from the same clusters, so 
there are similarities between all of these sequences.

[increasing hidden states](https://imgoat.com/uploads/6da2f590cd/176617.png)

Green, black and yellow  are "good","average" and "bad" sequences from the test set.
The model _does_ _not_ _know_ any of them. Further, they are dissimilar to train
set (beacuse of clustering done by CD-HIT). So we expect these _unknown_ _sequences_
to give poor scores. The average and bad sequenecs have poor scores, but 
surprisingly, with increased hidden states, the good sequences are scoring higher. 

This means, the model is detecting good sequences from entirely different and 
unknown sequences as well. (Not sure if this result is specific to our dataset.)

This means we can expect the model to generate *some* very good sequences (green),
even if it, for some reason, decides to give just the bad sequences (blue). 


# Generate
This generates sequences from given sequence. You have to specify the sequence,
model and number of sequence to generate. -v option recommended.


Example use case:
```console
$ python3 seqgen.py -s 
agtgtttgtgtctgcaatcccaagtttgtttgcgctgaaatatgcgatgctcaatgttatgatctgcgtactaagccgcagatcatagtgggaact 
-m example_model.pkl -n 200 -v

Generating sequences..
[=========================================================] 100% Time:  0:00:00
Exporting scores to csv..
```
