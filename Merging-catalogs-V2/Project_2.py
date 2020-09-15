#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 00:37:42 2020

@author: atilapaes

Project 2

Objective: incorporate the results from manual QC to the event catalog

"""
import pandas

#%% Data importing and preparation

# Import catalog of single events (no duplicate present)
catalog_single = pandas.read_csv('result_catalog_merged_single.csv',index_col='Unnamed: 0')
catalog_single['datetime']=pandas.to_datetime(catalog_single['datetime'])

# Import catalog for manual QC
catalog_manual = pandas.read_csv('result_catalog_merged_manual.csv',index_col='Unnamed: 0')
catalog_manual['datetime']=pandas.to_datetime(catalog_manual['datetime']) #Convert the strings generated before in the datetime format used by Pandas


#%% A fast method to demostrate that there is at least one MFA event in each block
set_blocks=set(catalog_manual.double_block.values)
set_MFA=set(catalog_manual[catalog_manual['source']=='MFA'].double_block.values)

if set_blocks==set_MFA:
    print('MFA in all blocks')

#%% With that in mind, let's keep just the MFA entries
catalog_manual_filtered=catalog_manual[catalog_manual['source']=='MFA'].copy()
catalog_manual_filtered['source']='Both'

catalog_manual_filtered=catalog_manual_filtered.drop(columns=['double_flag','manual','double_block'])
catalog_manual_filtered.reset_index(drop=True,inplace=True)

#%% Adding info from manual QC (a new event identified)

# Preparing a dataframe with the info from the new event
new_event = pandas.DataFrame({'datetime': ['2016-11-03 12:36:05'], 'Date': ['2016-11-03'], 'source':['ES']})
new_event['datetime']=pandas.to_datetime(new_event['datetime'])

# Merging the new event to the manual QC catalog
catalog_manual_filtered=pandas.concat([catalog_manual_filtered,new_event],ignore_index=True) #Concatenate the catalogs and ignore the index

#%% Merge single and manual filtered catalos
catalog_merged=pandas.concat([catalog_manual_filtered,catalog_single],ignore_index=True)

catalog_merged=catalog_merged.sort_values(by='datetime') # Sort values based in date and time
catalog_merged.reset_index(drop=True,inplace=True) # Reset index after sorting

#%%
catalog_merged.to_csv('result_catalog_merged_ES_MFA_FINAL.csv',sep=',',line_terminator='\n',index=True)