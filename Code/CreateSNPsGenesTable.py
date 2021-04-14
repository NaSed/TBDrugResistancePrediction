#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 20:31:10 2020

@author: nafiseh
"""

# This script creates a table of SNPs with position as index and having REF, ALT, Gene, Alias and DistanceFromBeginingOfGene columns.
import numpy as np
import pandas as pd
import pickle


#Home = '/Users/nafiseh/Google Drive/DrugResistance/'
Home = '/project/6004228/sedaghat/DR_Project/'
Data_path = Home + 'Data/'


##############################################
#           Creating SNP lists
##############################################

SNPList = pd.read_csv(Data_path + '/SNPs.csv', header = 0)
SNPList['POS'] = SNPList['POS'].astype(int)
type(SNPList.loc[8,'POS'])
##############################################
#           Reading reference genome
##############################################

refgene_df = pd.read_csv(Data_path + '/EPFL_Data/Mycobacterium_tuberculosis_H37Rv_allGenes.txt', sep='\t', header = 0, 
                    usecols=['Start', 'Stop', 'Strand', 'Locus', 'Name'], index_col=['Locus', 'Name'])

refgene_df = refgene_df.dropna(subset=['Start', 'Stop'])    
refgene_df['Start'] = refgene_df['Start'].astype(int)
refgene_df['Stop'] = refgene_df['Stop'].astype(int)

 
refgeneNames = refgene_df.index.values.tolist() # name


##############################################
#           Reading SNPs
##############################################

snp_df = SNPList.copy()

def FindSNPSOccurredInGenes(gene_start, gene_end, snp_lst):
    pos_arr = np.unique(list(snp_df['POS'].values))
    pos = list(pos_arr[np.where((gene_start < pos_arr) & (gene_end > pos_arr))])
    pos = sorted(pos, reverse = True)
    return pos

snp_df.insert(3, 'Locus', None)
snp_df.insert(4, 'Name', None)
snp_df.insert(5, 'DistanceFromBeginingOfGene', None)


for i in range(len(refgene_df)): # for each gene
    gene = refgene_df.index.values[i]
    strt, end, orientation = refgene_df.iloc[i]
    positions = FindSNPSOccurredInGenes(gene_start=strt, gene_end=end, snp_lst=snp_df)
#    snp_df['POS'].isin(positions)
    snp_df.loc[snp_df['POS'].isin(positions), 'Locus'] = refgeneNames[i][0]
    snp_df.loc[snp_df['POS'].isin(positions), 'Name'] = refgeneNames[i][1]
    for j in positions:
        snp_df.loc[snp_df['POS']==j, 'DistanceFromBeginingOfGene'] = j - strt

snp_df.to_csv(Data_path + 'SNPsGenes.csv', index=False)
pickle.dump(snp_df, open(Data_path + 'SNPsGenes.pkl', "wb"), protocol=pickle.HIGHEST_PROTOCOL)
