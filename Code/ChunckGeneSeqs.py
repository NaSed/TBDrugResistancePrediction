#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:49:46 2020

@author: nafiseh
"""
import pandas as pd
import pickle

#home_path = '/project/6004228/sedaghat/DR_Project/'
home_path = '/Users/nafiseh/Google Drive/DrugResistance/'
Data_path = home_path + 'Data/'
Result_path = home_path + 'Results/'

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Reading gene sequences
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def LoadSeqData(ind1, ind2, path):
    lst = list(range(ind1, ind2+1, 200))
    genes = list()
    isolates = list()

    for i in range(len(lst)-1):
        fileName = 'Isolate_Gene_Seqs_' + str(ind1)+ '_' + str(ind2) + '.pkl'
        print(fileName)

        dt = pickle.load(open(path + 'GeneSeqFiles/'+ fileName, "rb"))
        genes.append(dt.columns.to_list())
        isolates = dt.index.to_list()
        if i==0:
            geneSeqs = dt.copy()
        else:
            geneSeqs = pd.concat([geneSeqs, dt], axis='columns', ignore_index=True)
    print(geneSeqs.shape)
    print(len(genes))
    print(len(isolates))
    geneSeqs.index = isolates
    geneSeqs.columns = [item for sublist in genes for item in sublist]
    return(geneSeqs)


def expandGeneSeqFiles(dt, desired_length):
  # dt: a dataframe whose rows and columns correspond to isolates and genes, respectivly. Each entry is gene sequence.
  # desired_length: the desired length for each gene sequence
  # Note that if a gene sequence is longer than 'desired_length' that results in breaking the sequence to several parts (columns),
  # the name of columns would be (locus$i, gene$i) that i indicates to the index of part.

  def chunkstring(string, length):
    return [string[0+i:length+i] for i in range(0, len(string), length)]

  isolates = dt.index.to_list()
  genes = dt.columns.to_list()
  n1 = max([max(dt[genes[i]].apply(len)) for i in range(dt.shape[1])]) # gene length needed before
  # separating n-length chunks with '_' in each sequence
  for i in range(dt.shape[1]): # over genes
    dt[genes[i]] = dt[genes[i]].apply(lambda x: "_".join(chunkstring(x, desired_length)))


  newGeneSeqs = pd.DataFrame()
  for i in range(dt.shape[1]): # over genes
      res = dt[genes[i]].str.split('_', expand=True) # expanding a column to multiple columns

      if res.shape[1] > 1: # The max length string in the column is greater than the threshold (we have '_' in the string)
        newCol_n = len(list(res.columns)) # Number of new columns
        locus = dt.columns[i][0]
        name = dt.columns[i][1]

        res.columns =  [(locus + '$' + str(j), name + '$' + str(j)) for j in range(newCol_n)]
        res.index = isolates

        if i==0:
          newGeneSeqs = res
        else:          
            newGeneSeqs = pd.concat([newGeneSeqs, res], axis=1)
      else: # There would not be any changes to the columns since the gene sequences are shorter than the threshold
        res.columns = [genes[i]]  
        if i==0:
            newGeneSeqs = res        
        else:
            newGeneSeqs = pd.concat([newGeneSeqs, res], axis=1)
  
  # Removing chunks that have the same value across all isolates
  numofUniqueVals = list(newGeneSeqs.nunique())  # Number of unique values in each column
  rem_inds = [i for i in range(newGeneSeqs.shape[1]) if numofUniqueVals[i] ==1 ] # Columns that must be removed
  print('number of columns with only one unique value:', len(rem_inds))
  newGeneSeqs.drop(newGeneSeqs.columns[rem_inds], inplace=True, axis='columns')
  
  print('Total gene length BEFORE splitting:', n1 * len(genes))
  print('Total gene length AFTER splitting:', newGeneSeqs.shape[1] * desired_length)
  return(newGeneSeqs)
  
###################################################################
###################################################################
###################################################################
  
stp = 200
desired_length = 200
lst = list(range(0, stp+1, 200))

for i in range(len(lst)-1):
    geneSeqs = LoadSeqData(lst[i], lst[i+1], path=Data_path)
    chunckedSeqs = expandGeneSeqFiles(geneSeqs, desired_length)

    numbers_train = []
    for isolate in chunckedSeqs.index.to_list():
        tmp = getIntegerEncoded_acrossIsolates(chunckedSeqs, isolate, desired_length)
        numbers_train.append(tmp)

    numbers_train = np.array(numbers_train, dtype='int8')
    fileName = 'Chuncked' + str(desired_length) +'_Integer_Isolate_Gene_Seqs_' + str(lst[i])+ '_' + str(lst[i+1]) + '.pkl'
    
    pickle.dump(numbers_train, open(Data_path + 'GeneSeqFiles/'+ fileName, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
