import numpy as np
import progressbar
import itertools
import time
from itertools import accumulate


def helicopter(kwrg):
    if kwrg=='start':
        for c in itertools.cycle('/-\|'):
            print(c, end = '\r')
            time.sleep(0.2)

    else:
       pass


def codon_splitter(sequences,n, verbose = False):
    """Splits a list of sequences into a list of a list characters at 
    groups of n. 

    Args:
        sequences: List of sequence
        n: Group of chars to split from the sequence

    Returns:
        Nested list of splitted sequence for each sequence in given list

    Raises:
        ValueError if a your supplied list has just one sequence.
    
    
    Description:
        Sequences is a list of two or more individual sequences,and n may be
        set to 3. So function creates a list of 3 characters
        (codons)for each sequence. So the returned codons_list is a nested 
        list of two lists each list being a list of codons.
        
        
    """
    pbar = progressbar.NullBar(min_value=0, max_value=None)
    if verbose == True:
        print('Separating codons..')
        pbar = progressbar.ProgressBar(widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ',
                                            progressbar.ETA()], maxval=len(sequences))
    
    if len(sequences)<2:
        return ValueError("Too few sequences to train!")
        sys.exit()
    

    
    codons_list = []
    split_sequence = lambda x, n: [x[i:i+n] for i in range(0, len(x), n)]
    pbar.start()
    for i in range(len(sequences)):
        codons_list.append(split_sequence(sequences[i],n))
        pbar.update(i)
    pbar.finish()
    progressbar.streams.flush() 
    
    return codons_list

def create_mapping(codons_list,*argv):
    """Creates a mapping for unique items in codon_list
    
    Args:
        codon_list: List of codons
        *agrv: either 'encode' or 'decode'. (Optional)
        
    Returns:
        Dictionary with numerical encodings for unique items in supplied
        list
        
    Raises:
    ValueError: If *argv is not 'encode' or 'decode'
    
    Description:
         
    """
    
    flattned_list = [item for sublist in codons_list for item in sublist]
    codons_set = sorted(list(set(flattned_list)))
    
    codon_to_n = {codon:n for n,codon in enumerate(codons_set)}
    n_to_codon = {n:codon for n,codon in enumerate(codons_set)}
    
    for val in argv:
        if val not in ('encode','decode'):
            raise ValueError("agrs must be either 'encode' or 'decode' only.")
            break
        if val == 'encode':
            return codon_to_n
            break
        if val == 'decode':
            return n_to_codon
            break
    return codon_to_n, n_to_codon


def codons_to_num(codons_list,codon_to_n):
    """encodes codons from a given list in numerals
    
    Args:
        codon_list: List of codons
        
    Returns:
        
        
    Raises:
    
    Description:
         
    """ 
    
    #calculate depth of passed list
    depth = lambda L: isinstance(L, list) and (max(map(depth, L)) + 1) if L else 1
    encoded_list = []
    
    if depth(codons_list) == 1:
        encoded_list.append([codon_to_n[codon] for codon in codons_list])
    else:
        for i in range(len(codons_list)):
            encoded_list.append([codon_to_n[codon] for codon in codons_list[i]])

    return encoded_list


def num_to_codon(encoded_list,n_to_codon):
    """Decodes numerical values from a given list to codons
    
    Args:
        encoded_list: List of numerical encoded codons
        num_to_codons: A dictionary to decode numbers to codons. This can be output
                        of create_mapping with decode argv
        
    Returns:
        decoded list of numerals to codon sequence
        
    Raises:
    
    Description:
         
    """ 
    
    #calculate depth of passed list
    depth = lambda L: isinstance(L, list) and (max(map(depth, L)) + 1) if L else 1
    
    if depth(encoded_list) == 1:
            decoded_list.append([n_to_codon[num] for num in encoded_list])
    else:
        for i in range(len(encoded_list)):
            decoded_list.append([n_to_codon[num] for num in encoded_list[i]])
    return decoded_list


def convert_to_hmm_data(codons_list,codon_to_n,verbose = False):
    """Converts codons list to long list of sequences of observations (X).
    Before passing this to HMM make sure you pass the lengths of each sequence
    as a list, otherwise model thinks its a one big sequence, rather than a number
    of smaller sequences
    
    Args:
        codons_list: List of codons
        codon_to_num: A dictionary to encode codons to numbers. This can be output
                        of create_mapping with encode argv
        
    Returns:
        2D numpy-arrays of observations 
        
    Raises:
    
    Description:
         
    """ 
    pbar = progressbar.NullBar(min_value=0, max_value=None)
    if verbose == True:
        print('\nInitializing a model.')
        pbar = progressbar.ProgressBar(widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ',
                                            progressbar.ETA()], maxval=len(codons_list))

        
    #calculate depth of passed list
    depth = lambda L: isinstance(L, list) and (max(map(depth, L)) + 1) if L else 1
    
    

    X = []
  
    if depth(codons_list) == 1:
        num_list = []
        num_list.append([codon_to_n[codon] for codon in codons_list])
        X= np.array(num_list).T
    else:
        pbar.start()  
        for i in range(len(codons_list)):
            X.append([codon_to_n[codon] for codon in codons_list[i]])
            pbar.update(i)
        pbar.finish()
        progressbar.streams.flush() 
        
    return np.array(X).reshape(-1, 1)


def calc_scores(model,codons_list_for_hmm,lengths,verbose=False):
    
    pbar = progressbar.NullBar(min_value=0, max_value=None)
    if verbose == True:
        print('\nCalculating Scores..')
        pbar = progressbar.ProgressBar(widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                            progressbar.Percentage(), ' ',
                                            progressbar.ETA()], maxval=len(lengths)) 
       

        
    
    cumulative=list(accumulate(lengths))
    cumulative.insert(0,0) #to fix first range 0-32
    #scores = [mode.score(np.array(X[cumulative[length]:cumulative[length+1]]))]
    pbar.start()
    scores=[]
    for i in range(len(cumulative)-1):
        scores.append(model.score(codons_list_for_hmm[cumulative[i]:cumulative[i+1]]))
        pbar.update(i)
        
    pbar.finish()
    progressbar.streams.flush() 

    
    return scores



