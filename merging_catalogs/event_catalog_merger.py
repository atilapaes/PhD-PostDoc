#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 17:25:41 2019

@author: atilapaes


Merge the SRL and ES catalogs in a single filtered catalog of event date time
"""

import pandas, numpy


# Catalog srl and treatment
catalog_srl_initial=pandas.read_csv('catalog_srl.csv') # Read file
catalog_srl_initial['datetime']=catalog_srl_initial['Date']+'T'+catalog_srl_initial['Time'] # Conpose a datetime format readable by pandas
catalog_srl_initial['datetime']=pandas.to_datetime(catalog_srl_initial['datetime']) #Convert the strings generated before in the datetime format used by Pandas
catalog_srl_initial.drop(catalog_srl_initial.index[2023:],inplace=True) #Used to process just the first part of the dataset (first round of injection)
catalog_srl_initial['source']='srl' # Creating a column to be label that will identify the sorce identification method
    
# Catalog ES and treatment
catalog_es_initial=pandas.read_csv('catalog_es.csv') # Read file
catalog_es_initial['datetime']=pandas.to_datetime(catalog_es_initial['datetime']) # Convert column with string of datetime to the format of Pandas datetime
catalog_es_initial['date']=catalog_es_initial['datetime'].dt.date # Get a column with date only
catalog_es_initial['source']='es' # Creating a column to be label that will identify the sorce identification method

#%% Treatment of ES's catalog (there was some ducplicates)
catalog_es=catalog_es_initial[['datetime','source']].copy()
catalog_es.drop_duplicates(subset='datetime',inplace=True)
catalog_es.sort_values(by='datetime')#.reset_index()

#%% Treatment of SRL's catalog
catalog_srl=catalog_srl_initial[['datetime','source']].copy()

#%% Merge, sort and reset index
catalog=pandas.concat([catalog_es,catalog_srl],ignore_index=True) #Concatenate the catalogs and ignore the index
catalog=catalog.sort_values(by='datetime') # Sort values based in date and time

catalog.reset_index(drop=True,inplace=True) # Reset index after sorting

catalog['double']=False #Label for events identified by both methods

# The criterium for considering two consecutives events as the same is that they delay is =< 1.5 seconds
for index in range(1,len(catalog)):
    if (catalog.datetime[index] - catalog.datetime[index-1]) < pandas.Timedelta('1.5 second'):
        catalog.at[index-1,'source'] = 'Both' # Labelingone of the duplicated as "identified by both algorithms'
        catalog.at[index,'double'] = True # Labeling the second of duplicates as duplicate (to be removed following)

catalog=catalog[catalog.double ==False] # Droping (1 of 2) duplicated events 
catalog.reset_index(drop=True,inplace=True) #Reseting index for more readable catalog

catalog=catalog[['datetime','source']] # Exporting just the datetime and source catalog of each event


#%% Exporting the merged catalog
catalog.to_csv('catalog_merged.csv',sep=',',line_terminator='\n',index=False)
