import os
import subprocess
from subprocess import run, PIPE, Popen
import pandas as pd
import numpy as np



def access_calc(seq):
    os.environ["LD_PRELOAD"] = os.getcwd()+'/bin/loader.so'
    proc = run(['RNAplfold', '-W 210', '-u 210', '-O'], \
               stdout=PIPE,stderr=subprocess.DEVNULL,input=seq,\
               encoding='utf-8',env=os.environ) 
    return str(proc.stdout)
    
test_seq = 'ATGACGTCCAAAGTTTATGATCCGGAACAGCGCAAGAGGAATTAAAT'
_stdout = access_calc(test_seq)
openen = _stdout.split('%%EOF') #this is the openen values
##to do
##parsing the table
