#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:41:57 2020

@author: nafiseh
"""
import numpy as np
import pandas as pd

MappingTable = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }

###############################################################################
###############################################################################

def translate_TB_Gene(seq): 
# This function gets a sequence of 'A', 'C', 'G' and 'T' and translates it to a protein sequence

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582765/
# Start codons: GTG, ATG, TTG
# Stop codons: TAG, TAA, and TGA

# If start codons appear in the beginning they are translated to 'M',
# If they appear in the middle of sequence they translated to the corresponding amino based on the table

    
# Unit test: (reference: http://genome.tbdb.org/cgi-bin/GeneDetails.html?id=SRv3285)
#Rv3285 = 'GTG GCT AGT CAC GCC GGC TCG AGG ATC GCT (10) CGG ATC TCT AAG GTT CTC GTC GCC AAT CGC (20) GGC GAG ATC GCA GTG CGG GTG ATC CGG GCG (30) GCC CGC GAC GCC GGC CTG CCC AGC GTG GCG (40) GTG TAC GCCGAACCCGACGCCGAGTCCCCGCATGTTCGGCTGGCCGACGAGGCGTTCGCGCTGGGCGGCCAGACCTCGGCGGAGTCCTATCTGGACTTCGCCAAGATCCTCGACGCGGCAGCCAAGTCCGGGGCCAACGCCATCCACCCCGGCTACGGCTTCCTAGCGGAAAATGCCGACTTCGCCCAGGCGGTGATCGACGCCGGCCTGATCTGGATCGGCCCCAGCCCGCAGTCGATCCGCGACCTGGGCGACAAGGTCACGGCCCGTCACATCGCGGCCCGCGCTCAGGCGCCCCTGGTGCCGGGTACCCCCGATCCGGTCAAAGGCGCCGACGAGGTGGTGGCATTCGCCGAGGAGTACGGCCTGCCGATCGCGATCAAGGCCGCCCACGGCGGCGGCGGCAAGGGCATGAAGGTGGCCCGCACCATCGACGAGATTCCGGAGCTGTACGAGTCGGCGGTGCGCGAGGCCACGGCCGCGTTCGGCCGCGGTGAGTGCTACGTGGAGCGCTATCTCGACAAGCCGCGCCACGTCGAAGCACAGGTGATCGCCGACCAGCACGGCAACGTCGTCGTCGCCGGCACCCGGGACTGCTCGCTGCAGCGCCGCTACCAGAAGCTGGTCGAGGAGGCGCCCGCACCGTTCCTGACCGACTTTCAACGCAAAGAGATCCACGACTCGGCCAAACGGATTTGCAAAGAGGCCCATTACCACGGCGCCGGCACCGTCGAATACCTGGTCGGTCAGGACGGCTTGATCTCGTTCTTGGAGGTCAACACGCGCCTTCAGGTAGAACACCCGGTCACCGAGGAAACCGCGGGCATCGACTTGGTGCTGCAGCAATTCCGGATCGCCAACGGCGAAAAGCTGGACATCACCGAGGATCCCACCCCGCGCGGGCACGCCATCGAATTCCGGATCAACGGCGAGGACGCGGGGCGTAACTTCCTACCGGCGCCCGGGCCGGTGACAAAGTTCCACCCGCCGTCCGGCCCCGGTGTGCGGGTGGACTCCGGTGTCGAGACCGGCTCGGTGATCGGCGGCCAGTTCGACTCGATGCTGGCCAAGCTGATCGTGCACGGTGCCGACCGCGCCGAGGCGCTGGCGCGGGCCCGGCGCGCGCTGAACGAGTTCGGTGTCGAAGGCCTGGCGACGGTCATCCCGTTTCACCGCGCCGTGGTGTCCGACCCGGCATTCATCGGCGACGCGAACGGCTTTTCGGTACATACCCGCTGGATCGAGACCGAGTGGAATAACACCATCGAGCCCTTTACCGACGGCGAACCTCTCGACGAGGACGCCCGGCCGCGTCAGAAGGTGGTCGTCGAAATCGACGGTCGCCGCGTCGAAGTCTCGCTGCCGGCTGATCTCGCGCTGTCCAATGGCGGCGGTTGCGACCCGGTCGGTGTCATCCGGCGCAAGCCCAAGCCGCGCAAGCGGGGTGCGCACACCGGCGCGGCGGCCTCCGGTGACGCGGTGACCGCGCCTATGCAGGGCACCGTAGTTAAGTTCGCGGTCGAAGAAGGGCAAGAGGTCGTGGCCGGCGACCTAGTGGTG'+'GTCCTCGAGGCGATGAAGATGGAAAACCCGGTCACCGCGCATAAGGATGGCACCATCACC'+'GGGCTGGCGGTCGAGGCGGGCGCGGCCATCACCCAGGGCACGGTGCTCGCCGAGATCAAGTAA'
#prSeq, stop, premature = translate_TB_Gene(Rv3285)
#true_prSeq = 'MASHAGSRIA (10) RISKVLVANR (20) GEIAVRVIRA (30) ARDAGLPSVA (40) VYAEPDAESPHVRLADEAFALGGQTSAESYLDFAKILDAAAKSGANAIHPGYGFLAENADFAQAVIDAGLIWIGPSPQSIRDLGDKVTARHIAARAQAPLVPGTPDPVKGADEVVAFAEEYGLPIAIKAAHGGGGKGMKVARTIDEIPELYESAVREATAAFGRGECYVERYLDKPRHVEAQVIADQHGNVVVAGTRDCSLQRRYQKLVEEAPAPFLTDFQRKEIHDSAKRICKEAHYHGAGTVEYLVGQDGLISFLEVNTRLQVEHPVTEETAGIDLVLQQFRIANGEKLDITEDPTPRGHAIEFRINGEDAGRNFLPAPGPVTKFHPPSGPGVRVDSGVETGSVIGGQFDSMLAKLIVHGADRAEALARARRALNEFGVEGLATVIPFHRAVVSDPAFIGDANGFSVHTRWIETEWNNTIEPFTDGEPLDEDARPRQKVVVEIDGRRVEVSLPADLALSNGGGCDPVGVIRRKPKPRKRGAHTGAAASGDAVTAPMQGTVVKFAVEEGQEVVAGDLVVVLEAMKMENPVTAHKDGTITGLAVEAGAAITQGTVLAEIK_'
#prSeq == true_prSeq


    if seq == np.NaN:
        return None
     
    protein ="M" 
    
    startInd = [seq.find('GTG'), seq.find('TTG'), seq.find('ATG')]
    # if start codon does not exist find(.) returns -1
    if len(np.unique(startInd))==1 and np.unique(startInd) == -1:
        return None
    startInd = min([x for x in startInd if x !=-1 ])

    assert(type(startInd) == np.int)
    prematureStop = False
    for i in range(startInd+3, len(seq), 3): 
        if i+3 >= len(seq):
            print('The remaining part of sequence is not translated!')
            break
        codon = seq[i:i + 3]
        
        if MappingTable[codon] == '_':
            if (i+3) < (len(seq)-1):
                prematureStop = True
            break
        protein += MappingTable[codon]

    stopInd = i    

    return [protein, stopInd, prematureStop]

###############################################################################
###############################################################################
        
def sequence_complement(seq):
    # It gets a sequence of 'A', 'C', 'G' and 'T' and returns its complement
    compSeq = []
    for i in range(len(seq)):
        if seq[i] == 'T':
            compSeq.append('A')
        if seq[i] == 'A':
            compSeq.append('T')
        if seq[i] == 'C':
            compSeq.append('G')
        if seq[i] == 'G':
            compSeq.append('C')
    return(compSeq)

###############################################################################
###############################################################################
        
       
def gene_sequence(strt, end, orientation, refgenome):
    # It gets start and end positions as well as orientation of a gene along with 
    # reference genome and returns sequence of the gene considering its orientation.
    if orientation == '+':
        return(refgenome[(strt-1):end])
    else:
        seq = refgenome[(strt-1):end]
        seq.reverse()
        seq = sequence_complement(seq)
        return(seq)
        
###############################################################################
###############################################################################
        
# This function is called after translation and returns premature flag
def getPrematureFlag(lst):
    if lst == None:
        return None
    else:
        return lst[2]
###############################################################################
###############################################################################
        
# This function is called after translation and returns protein sequence
def getProteinSeq(lst):
    if lst == None:
        return None
    else:
        return lst[0]
###############################################################################
###############################################################################
        
# This function is called after translation and returns position of stop codon
def getStopInd(lst):
    if lst == None:
        return None
    else:
        return lst[1]

###############################################################################
###############################################################################

# This function gets a sequence and applies SNPs on it.
# seq: it is a string
# snps: it is a dataframe whose index indicates to the position of SNP and it has two columns REF and ALT which 
# REF shows the letter in the original sequence and ALT shows the letter that REF is replaced with.
# It returns the sequence after changing the SNPs
def seqAfterSNPs(seq, snps):
    snps.sort_index(inplace=True, ascending=False)
    seq = list(seq)
    positions = snps.index.values.tolist()
    for i in range(len(snps)):
        REF, ALT = snps.loc[positions[i],['REF', 'ALT']]
        ind = positions[i]
        assert("".join(seq[ind:(ind+len(REF))]) == REF)
        a=seq[0:ind]
        a.extend(ALT)
        a.extend(seq[ind+len(REF):])
        seq = a

    return(seq)
    
###############################################################################
###############################################################################

# It gets a list of SNPs positions (Distance from the begining of the gene) and a gene sequence.
# It returns positions of SNPs whose occurrence are after a premature stop codon
def SNPsAfterStopCodon(snpsPos, geneSeq):
    translated = translate_TB_Gene(geneSeq)
    prematureFlag = translated[2]
    stopPos = translated[1]
    if prematureFlag == False:
        return (None)
    else:
        afterStop = [x for x in snpsPos if x > stopPos]
        return (afterStop)
    


###############################################################################
###############################################################################

# This function determines type of mutations: missense, nonsense, silent or frameshift
# and add a column 'Type' to 'snps' data frame.
#
#refgeneSeq: it is string, the sequence that snps are applied on
#snps: This is a data frame with columns DistanceFromBeginingOfGene, ALT and REF whose indices are positions of SNPs from the beginning of sequence

def FindMutationType(refgeneSeq, snps):
#    snps = SNPs.copy()
    refgeneSeq = list(refgeneSeq)
    assert(len(snps) > 0), "There is no SNPs"
    
    refCol = np.where(snps.columns=='REF')[0][0]
    distCol = np.where(snps.columns=='DistanceFromBeginingOfGene')[0][0]
    altCol = np.where(snps.columns=='ALT')[0][0]
    
    snps['change'] = None
    
    
    # Just for testing correctness
    for i in snps.index:
        REF, ALT = snps.loc[i, ['REF', 'ALT']]
        ind = snps.loc[i, 'DistanceFromBeginingOfGene']
        if "".join(refgeneSeq[ind:(ind+len(REF))]) != REF:
            assert("".join(refgeneSeq[ind:(ind+len(REF))]) == REF), 'Just for testing match between refgenome and REF values'
        snps.loc[i, 'change'] = len(ALT) - len(REF)

    indices = snps.index.tolist()

    changeVec = snps['change'].values.tolist()
    nonzero_ind = np.nonzero(changeVec)[0].tolist()
    
    tempSum = 0
    
    # Frameshift regions
    strt_ind = []
    end_ind = []
    for i in nonzero_ind:
        ch = changeVec[i]
        if tempSum == 0 and ch%3 != 0:
            strt_ind.append(i)
        if ch % 3 !=0:
            tempSum = tempSum + ch
            if tempSum == 0:
                end_ind.append(i)
    
    
    if len(end_ind) < len(strt_ind):
        end_ind.append(len(changeVec)-1)
    
    snps['Type'] = None
    for i in range(len(strt_ind)):
        if len(strt_ind)>1:
            snps.loc[indices[strt_ind[i]]:indices[end_ind[i]], 'Type'] = 'Frameshift'
    
    snp_types = snps.Type.values.tolist()
    if None not in snp_types:
        return(snps)
        
    
    for i in range(len(snps)):
        if snp_types[i] == None:
            ref = snps.iloc[i, refCol]
            alt = snps.iloc[i, altCol]
            
            if len(ref) < len(alt): # insertion
                if (len(alt) - len(ref))%3 != 0: # Three letters are inserted, so the other
                    snp_types[i] = 'Insertion'
                else:
                    snp_types[i] = 'Frameshift'
            elif len(ref) > len(alt): # deletion
                if (len(ref) - len(alt))%3 != 0: # Three letters are deleted, so the other
                    snp_types[i] = 'Deletion'
                else:
                    snp_types[i] = 'Frameshift'
            else: # One point mutation
                diffCount = sum(1 for a, b in zip(ref,alt) if a != b) # The number of different chars (incase both ref and alt have len>1)
                assert(diffCount < 2), 'it is identified as point mutation but there is more than one change in the ref'
                chk= False
                if len(ref) > 1 and len(alt)>1:
                    refl = list(ref)
                    altl = list(alt)
                    for kk in range(1,len(refl)):
                        if ref[-kk]==alt[-kk]:
                            print(kk)
                            del refl[-1]
                            del altl[-1]
                        else:
                            break
                    chk = False
                    if len(refl)>1 and len(altl)>1:
                        for kk in range(len(refl)):
                            if refl[kk]!= altl[kk]:
                                chk = True
                                break
                    if chk == True:
                        dist = snps.loc[indices[i], 'DistanceFromBeginingOfGene']-kk
                        refl = refl[kk:]
                        altl = altl[kk:]
                    ref = ''.join(refl)
                    alt = ''.join(altl)
                    
                pos = int(snps.iloc[i, distCol] % 3)
                if chk == True: # REF and ALT had len>1 while the first chars are the same
                    ind = dist
                else:
                    ind = snps.loc[indices[i], 'DistanceFromBeginingOfGene']
                if pos == 0: # The position that SNP has occurred is the first letter of a three-letter codon
                    before_change = refgeneSeq[ind:ind+3]
                    after_change = before_change.copy()
                    after_change[pos] = alt
                elif pos == 1: # The position that SNP has occurred is the second letter of a three-letter codon
                    before_change = refgeneSeq[(ind-1):(ind+2)]
                    after_change = before_change.copy()
                    after_change[pos] = alt
                else:  # The position that SNP has occurred is the third letter of a three-letter codon
                    before_change = refgeneSeq[(ind-2):(ind+1)]
                    after_change = before_change.copy()
                    after_change[pos] = alt
                   
                # Mapping the codons to amino acids    
                aft = MappingTable[''.join(after_change)]
                bef = MappingTable[''.join(before_change)]
                
                if aft == '_':
                    snp_types[i] = 'Nonsense'
                elif bef == aft:
                    snp_types[i] = 'Silent'
                elif bef != aft:
                    snp_types[i] = 'Missense'
    snps['Type'] = snp_types
    snps.drop(['change'], axis = 'columns', inplace = True)
    return snps

###############################################################################
###############################################################################

