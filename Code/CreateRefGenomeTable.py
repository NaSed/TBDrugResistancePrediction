#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 12:44:19 2020

@author: nafiseh
"""

# This scripts creates reference genome table with columns: 
# 'start', 'end', 'orientation', 'genSeq', 'proteinSeq', 'translationStopInd', 'prematureStop'

import numpy as np
import pandas as pd
import DNA_Manipulation.utils as DM
import pickle 

Home = '/project/6004228/sedaghat/DR_Project/Data/'

#Home = '/Users/nafiseh/Google Drive/DrugResistance/Data/'

# Files we already have
RefGenomeFile = Home  + 'EPFL_Data/Mycobacterium_tuberculosis_H37Rv_genes_v3.csv'     

# New file to save
RefGeneTableFile = Home + 'EPFL_Data/RefGeneTable.pkl'

########################################################
########################################################

geneSeqs = pd.read_csv(RefGenomeFile, header = 0, usecols=['Start', 'Stop', 'Strand', 'Locus', 'Name', 'Sequence', 'Type'], index_col=['Locus', 'Name'])
geneSeqs = geneSeqs.loc[geneSeqs.Type == 'CDS',]
geneSeqs.drop(['Type'], inplace = True, axis = 'columns')
refgeneNames = geneSeqs.index.values.tolist() # name
geneSeqs['OriginalSeq'] = None # The original sequence before orientation
for gene in refgeneNames:
    if geneSeqs.loc[gene,'Strand'] == '-': # If strand is minus we convert to the real reference genome sequence  .values[0]
        refseq = list(geneSeqs.Sequence.loc[gene,]) # It is real sequence of gene  .values[0]
        refseq.reverse()
        geneSeqs.loc[gene, 'OriginalSeq'] = ''.join(DM.sequence_complement(refseq))
    else:
        geneSeqs.loc[gene, 'OriginalSeq'] = geneSeqs.loc[gene, 'Sequence']
        

geneSeqs['ProteinSeq'] = None
geneSeqs['translationStopInd'] = None
geneSeqs['prematureStop'] = None

for gene in refgeneNames:
    res = DM.translate_TB_Gene(geneSeqs.loc[gene, 'Sequence'])
    geneSeqs.loc[gene, 'ProteinSeq'] = res[0]
    geneSeqs.loc[gene, 'translationStopInd'] = res[1]
    geneSeqs.loc[gene, 'prematureStop'] = res[2]

# Saving the file
pickle.dump(geneSeqs, open(RefGeneTableFile, "wb"), protocol=pickle.HIGHEST_PROTOCOL)

#ref_df.to_csv(Home + 'RefGenome/RefGeneTable.csv', index=True)
