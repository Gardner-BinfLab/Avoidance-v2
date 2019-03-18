# Optimization 

This script selects and possibly further optimizes the sequences generated from a trained Markov model. The details on how train a Markov model and generate synonymous sequences are [here.](https://github.com/Gardner-BinfLab/Avoidance-v2/tree/master/Markov_model)

## Working
![working](https://i.imgur.com/chjaaAO.png)

## Arguments
**required args**

```-m MRNA ``` : mRNA sequences in fasta or csv. If sequences are in csv, please ensure your sequences are in a column named 'sequence'.

```-r RANDOMFOREST ``` : trained random forest model. Details on how to do so are [here.](https://github.com/Gardner-BinfLab/Avoidance-v2/tree/master/Random_forest)

```[-u UTR5]``` : 5′ UTR (71 nucleotides only!). Deafult is the pET vector 'aggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacc'. 

**optional args**

```[-s]``` : this switches simulated annealing (SA) on. Default: no SA

```[-c COUNT]``` : number of top sequences to pick for further optimization by SA. Default 10.

```[-g GEN]``` : number of new sequences to generate per sequence. Default: 1

```[-n NITER]``` : number of iterations to do in a SA. Default: 200

```[-o OUTPUT]``` : output file name. default 'optimized\_sequences'. The full analysis of the given mRNA sequences and the chosen sequences are  also exported as 'output\_name\_mrna_analysis' and 'chosen' respectively to `results\optimized_sequences\` folder. 


**example usage :**

```console
bash-4.2$ python3 optimizer.py -m test_sequences.csv -r rf_expr.sav -o test_output -s
================================================
20190316-103303
using  aggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacc  as the 5' utr..
calculating features for given sequence list..
done!
loading random forest model..
done!
selecting 10 good sequences..
the max predicted probability of sequences being highly expressed is  0.955
generating 1000 random synonymous sequences for a background model..
calculating features for background sequences..
optimization started..this may take a while..
|██████████████████████████████████████████████████| 100% (2/2) at sequence :2
Completed!
optimization completed! exporting sequences...
full process was completed successfully!
20190316-103631
```

