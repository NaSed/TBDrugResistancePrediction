#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jun 12  3 15:52:00 2020

@author: nafiseh
"""

import DNA_Manipulation.utils as DM
import pickle 

Home = '/Users/nafiseh/Google Drive/DrugResistance/Data/'
#Home = '/project/6004228/sedaghat/DR_Project/Data/'

# Gene and protein sequences from reference genome created by CreateRefGenomeTable.py
RefGeneTableFile = Home + 'EPFL_Data/RefGeneTable.pkl'

strt = 0
stp = 4200
########################################################
#
#         Reading genes and proteins in reference genome 
#
########################################################
ref_df = pickle.load(open(RefGeneTableFile, 'rb'))

# Removing genes whose sequences have premature stop codons
ref_df = ref_df[ref_df['prematureStop'] == False]
ref_df.shape

# Removing those genes whose protein sequences are very small
ref_df = ref_df[ref_df.ProteinSeq.apply(len) > 20]
    

inds = list(range(strt, stp+1, 200))

for i in range(0, len(inds)-1):

    geneSeqsFile = Home + 'GeneSeqFiles/Isolate_Gene_Seqs_' + str(inds[i]) + '_' + str(inds[i+1]) + '.pkl'
    
    # New files to save
    proteinSeqsFile = Home + 'ProteinSeqFiles/Isolate_Protein_Seqs_' + str(inds[i]) + '_' + str(inds[i+1]) + '.pkl'
    prematureFile = Home + 'PrematureFiles/Premature_Gene_' + str(inds[i]) + '_' + str(inds[i+1]) + '.pkl' 
    percentTranslatedFile = Home + 'percentTranslatedFiles/translationPercentage_Gene_' + str(inds[i]) + '_' + str(inds[i+1]) + '.pkl' 

    ######################################################################
    #
    #                           TB Isolates
    #
    ######################################################################
    geneSeqs = pickle.load(open(geneSeqsFile,'rb'))
    
    # Removing genes that are not from MTb and have premature stop codon in the reference genome
    mustKeep = set(ref_df.index.tolist())
    allGenes = set(geneSeqs.columns.tolist())
    mustKeep = list(allGenes.intersection(mustKeep))    
    geneSeqs = geneSeqs.filter(mustKeep, axis='columns')
    
    # Translating to protein sequences, finding premature stop codons, and translation stop index
    translated = geneSeqs.applymap(DM.translate_TB_Gene)
    
    premature = translated.applymap(DM.getPrematureFlag)
    proteinSeqs = translated.applymap(DM.getProteinSeq)
    stopInd = translated.applymap(DM.getStopInd)
    
    TranslatedPercent = stopInd.copy()
    seqlen = geneSeqs.applymap(len)
    TranslatedPercent = TranslatedPercent.div(seqlen)
    
    # Removing genes that have the same value across all isolates
    #    numofUniqueVals = list(stopInd.nunique())  # Number of unique values in each column
    #    rem_inds = [i for i in range(stopInd.shape[1]) if numofUniqueVals[i] ==1 ] # Columns that must be removed
    #    stopInd.drop(stopInd.columns[rem_inds], inplace=True, axis='columns')
    
    pickle.dump(premature, open(prematureFile, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(proteinSeqs, open(proteinSeqsFile, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(TranslatedPercent, open(percentTranslatedFile, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    
    print("DONE")
