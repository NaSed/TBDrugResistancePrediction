#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:49:59 2020

@author: nafiseh
"""

# This script concatenate small dataframes stored in a directory
import glob
import pandas as pd
import pickle

Home = '/Users/nafiseh/Google Drive/DrugResistance/Data/'
Data_path = Home+ 'translationPercentage/'
mylist = [f for f in glob.glob(Data_path + "*.pkl")]
mylist = sorted(mylist)
fileId = 'translationPercentage_Gene'
strt = 0
stp = 4200
lst = range(strt, stp+1, 200)
dt = pickle.load(open(Data_path+ fileId+'_'+str(lst[0])+'_'+str(lst[1])+'.pkl', 'rb'))
for i in range(1, len(lst)-1):
        name = Data_path + fileId + '_'+str(lst[i])+'_'+str(lst[i+1])+'.pkl'
        print(name)
        y = pickle.load(open(name, 'rb'))
        print(y.shape)
        dt = pd.concat([dt,y], ignore_index=False, axis = 'columns')

pickle.dump(dt, open(Data_path+ fileId+'_'+str(0)+'_'+str(3990)+'.pkl', 'wb'))


#dt = pd.concat([
#pd.concat([pd.Series(['Training'] * len(train_isolates)), pd.Series(train_isolates)], axis = 1),
#pd.concat([pd.Series(['Validation'] * len(val_isolates)), pd.Series(val_isolates)], axis = 1),
#pd.concat([pd.Series(['Test'] * len(test_isolates)), pd.Series(test_isolates)], axis = 1)], axis = 0)
#
#
#dt.columns = ['Type', 'Isolate']
#
#dt.to_csv(Data_path+'Training_validation_test_Data.csv',index=False)
