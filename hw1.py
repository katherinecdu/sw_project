#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/python
__author__ = "Katherine Du"
__email__ = "katherine.du@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
from numpy import *
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,
default=-2, type=int)
parser.add_argument('-e', '--extgap', help='extension gap', required=False,
default=-1, type=int)
args = parser.parse_args()

arr_seq = np.array([])
f1=open(args.input,'r')
for line in f1.readlines():
    arr_seq = np.append(arr_seq, line.split())
seq2 = arr_seq[0]
seq1 = arr_seq[1]
    
arr = np.array([])
f2=open(args.score,'r')
for line in f2.readlines():
    arr = np.append(arr, line.split())
    
arr = ['', *arr]
arr = np.reshape(arr, (-1, 24))
BLOSUM62=arr

def match_score(alpha,beta): #def match_score(alpha,beta,BLOSUM62):
    alphabet={}
    alphabet["A"] = 1
    alphabet["B"] = 2
    alphabet["C"] = 3
    alphabet["D"] = 4
    alphabet["E"] = 5
    alphabet["F"] = 6
    alphabet["G"] = 7
    alphabet["H"] = 8
    alphabet["I"] = 9
    alphabet["K"] = 10
    alphabet["L"] = 11
    alphabet["M"] = 12
    alphabet["N"] = 13
    alphabet["P"] = 14
    alphabet["Q"] = 15
    alphabet["R"] = 16
    alphabet["S"] = 17
    alphabet["T"] = 18
    alphabet["V"] = 19
    alphabet["W"] = 20
    alphabet["X"] = 21
    alphabet["Y"] = 22
    alphabet["Z"] = 23
   
    lut_x=alphabet[alpha]
    lut_y=alphabet[beta]
    return int(BLOSUM62[lut_x][lut_y])

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    
    m,n =  len(seq1),len(seq2) #length of two sequences

    penalty=0;   #define the gap penalty
    opengap=openGap
    extgap=extGap

    #generate DP table and traceback path pointer matrix
    score=zeros((m+1,n+1))   #the DP table
    pointer=np.zeros([m+1,n+1,2])  #to store the traceback path
    
    max_score=0;
    all_gap_up = []
    all_gap_down = []
    score_up=0
    score_down=0

    #calculate DP table and mark pointers
    for i in range(1,m+1):
        for j in range(1,n+1):
        
            # gap up
            for val in range(1, i+1):
                all_gap_up.append(score[val][j] + opengap + (i-val-1)*extgap)
            score_up = max(all_gap_up)
            index_up = np.argmax(all_gap_up)
            #index_up = len(all_gap_up) - all_gap_up[::-1].index(score_up) # find index of last occurence of max val + 1
            all_gap_up=[]
        
            # gap down
            for val in range(1, j+1):
                all_gap_down.append(score[i][val] + opengap + (j-val-1)*extgap)
            score_down = max(all_gap_down)
            index_down = np.argmax(all_gap_down)
            #index_down = len(all_gap_down) - all_gap_down[::-1].index(score_down) # find index of last occurence of max val + 1
            all_gap_down=[]

            score_diagonal=score[i-1][j-1]+match_score(seq1[i-1],seq2[j-1]);
            score[i][j]=max(0,score_up,score_down,score_diagonal);

            if score[i][j]==0:
                pointer[i][j]=[0,0]; #0 means end of the path
            elif score[i][j]==score_diagonal:
                pointer[i][j]=[i-1,j-1]; #3 means trace diagonal
            elif score[i][j]==score_up:
                pointer[i][j]=[index_up+1,j] #1 means trace up
            elif score[i][j]==score_down:
                pointer[i][j]=[i,index_down+1] #2 means trace left
            if score[i][j]>=max_score:
                max_i=i;
                max_j=j;
                max_score=score[i][j];
       
    #END of DP table
    
    #alignments
    
    align1,align2='',''; #initial sequences

    i,j=max_i,max_j; #indices of path starting point

    pointer = pointer.astype(int)

    #traceback, follow pointers
    while i!=0 or j!=0:
        if i==0 or j==0:
            break
        if i-1==pointer[i,j][0] and j-1==pointer[i,j][1]:
            align1=align1+seq1[i-1];
            align2=align2+seq2[j-1];
        elif i==pointer[i,j][0]:
            for val in range(1,j-pointer[i,j][1]+1):
                align1=align1+'-';
                align2=align2+seq2[j-val];
        elif j==pointer[i,j][1]: #trace up, col same
            for val in range(1,i-pointer[i,j][0]+1):
                align1=align1+seq1[i-val];
                align2=align2+'-';
        if pointer[i,j][0] != 0 or pointer[i,j][1] != 0:
            endpoint = pointer[i,j]
        new_i=pointer[i,j][0]
        new_j=pointer[i,j][1]
        i=new_i
        j=new_j
    
    #END of traceback

    align1=align1[::-1]; #reverse sequence 1
    align2=align2[::-1]; #reverse sequence 2
    
    if endpoint[0]==1 and endpoint[1]==1:
        align2 = '(' + align2 + ')'
        align1 = '(' + align1 + ')'
    else:
        align2 = seq2[0:endpoint[1]] + '(' + align2 + ')' + seq2[max_j:len(seq2)]
        align1 = seq1[0:endpoint[0]] + '(' + align1 + ')' + seq1[max_i:len(seq1)]

    long = max(len(align2),len(align1))

    if len(align2) == long:
        align1 = " " * (endpoint[1]-endpoint[0]) + align1
        align1 = align1 + " " * (long-len(align1))
    elif len(align1) == long:
        align2 = " " * (endpoint[0]-endpoint[1]) + align2
        align2 = align2 + " " * (long-len(align2))
    
    mid = ''
    for val in range(0,long):
        if align2[val] == align1[val] and align2[val] != '(' and align2[val] != ')':
            mid = mid + '|'
        else:
            mid = mid + ' '
        
    # prepare table for output
    
    arr_seq1 = np.array([''])
    arr_seq2 = np.array(['',''])

    for val in seq2:
        arr_seq1 = np.append(arr_seq1, val)
    
    for val in seq1:
        arr_seq2 = np.append(arr_seq2, val)

    a = np.array(score.astype(int), dtype=object)
    a = np.insert(a, 0, arr_seq1, 0)
    a = np.insert(a, 0, arr_seq2, 1)
    
    ### write output
    out = open("output.txt", "w")

    out.write('-----------'+'\n')
    out.write('|Sequences|'+'\n')
    out.write('-----------'+'\n')
    out.write('sequence1'+'\n')
    out.write(seq2+'\n')
    out.write('sequence2'+'\n')
    out.write(seq1+'\n')
    out.write('--------------'+'\n')
    out.write('|Score Matrix|'+'\n')
    out.write('--------------'+'\n')

    for i in range(0,len(arr_seq2)):
        for j in range(0,len(arr_seq1)+1):
            out.write(str(a[i,j])+'\t')
        out.write('\n')

    out.write('----------------------'+'\n')
    out.write('|Best Local Alignment|'+'\n')
    out.write('----------------------'+'\n')
    out.write('Alignment Score:'+str(int(max_score))+'\n')
    out.write('Alignment Results:'+'\n')
    out.write(align2+'\n')
    out.write(mid+'\n')
    out.write(align1+'\n')

    out.close()
    
### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)

