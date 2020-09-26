#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 19:18:12 2020

@author: atilapaes


Merger of catalog of potential events

This script merges all CSVs in a specified folder, merge them, sort by date and export 
a CSAV file to be used as a 
"""
import pandas, glob


csv_folder='catalog_events/'
catalog_pe_folder='catalog_pe/'

# List the csv files in the folder
list_csv_pe=glob.glob(csv_folder+'*.csv')

# Read and concatenate the csv files
catalog_pe=pandas.read_csv(list_csv_pe[0])

for index_csv in range(1,len(list_csv_pe)):
    catalog_aux=pandas.read_csv(list_csv_pe[index_csv])
    catalog_pe=pandas.concat([catalog_pe,catalog_aux],ignore_index=True) #Concatenate the catalogs and ignore the index

#%% Catalog treatment

# Convert column to datetime
catalog_pe['datetime']=pandas.to_datetime(catalog_pe['datetime'])

# erase content in true_positve column
catalog_pe['true_positive']=''

# Sort by datetime
catalog_pe = catalog_pe.sort_values(by='datetime',ignore_index=True)

#%% Export catalog of potential events
catalog_pe.to_csv(catalog_pe_folder+'catalog_pe.csv',sep=',',line_terminator='\n',index=False)
