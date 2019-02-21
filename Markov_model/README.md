# Markov model
### Changelog
 - Major rewrite and refactoring. 
 - We now use the highly scalable Pandas dataframes instead of plain numpy arrays.
 - Useless dependencies removed.



# Install dependencies 
Requires ~`Python3.6+`~`Python3.4+`

This command installs the required dependencies.

`pip3 install --user scipy numpy pandas `


# Manual
  - `-h` gives possible options.
  - An example model is included and some training sequences are in test folder.


# Train
Train takes a list of sequences in .csv or .fasta and builds the model. The built model is 
exported as a python pickle (.pkl) file in 'results/model' subdirectory.
We train for 32 codons (96 nucleotides) by default but you can pass the length
in codons or just use flag `-a` to train over full length (not recommended).


```console
$ python3 train.py -f /directory/train_set_.csv 
Reading codons..
[===============                                   ] 30% (410/1333)
Stop codon encountered somewhere before the last position.
[==================================================] 100% (1333/1333)Finished!
1  sequence(s) were skipped because we encountered stop codons.

Training started. It may take a while..

Zeroth order
[==================================================] 100% (31/31)Finished!

Exporting model data to file..

```

# Scores
It calculates bitscores of seqence. You need to provide foreground and background model.
if a background is not provided, it will automatically make and use a uniform background.


Example use case:
```console
$ python3 score.py -f /dir/test.csv -m /dir/foreground.pkl -b /dir/backgnd.pkl
[===============                                   ] 30% (410/1334)Finished!
Stop codons encountered before the last position for sequence :  411
[================================================= ] 99% (1333/1334)Finished!
Exporting results..

```



# Generate
This generates sequences from given sequence. You have to specify the sequence,
model and number of sequences to generate.


Example use case:
```console
$ python3 seqgen.py -s 
agtgtttgtgtctgcaatcccaagtttgtttgcgctgaaatatgcgatgctcaatgttatgatctgcgtactaagccgcagatcatagtgggaact 
-m example_model.pkl -n 200 

Generating sequences..
[=========================================================] 100% (200/200)Finished!
Exporting scores to csv..
```
