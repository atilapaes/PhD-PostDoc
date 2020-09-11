#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 18:27:50 2020

@author: atilapaes

Merging catalogs V2
This script is intended to provide a more reliable method to merge two different catalogs.
Also, the extra information present in SRL dataset is preserved


Part 1 - Similar to the previously developed Proj01_catalog_merger_double_with_MW.py)

Part 2 - Split the merged catalog into single and doubles

Part 3 - For the double catalog: plot the waveforms and time for the events, export the figures for manual inspection

"""

import pandas, obspy, numpy
from tqdm import tqdm
from obspy.core import UTCDateTime

#%%### Part 1 #################################################################
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

# Treatment of ES's catalog (there was some ducplicates)
catalog_es=catalog_es_initial[['datetime','source']].copy()
catalog_es.drop_duplicates(subset='datetime',inplace=True)
catalog_es.sort_values(by='datetime')#.reset_index()

# Treatment of SRL's catalog. The 3D coordinates and the MW are preserved
catalog_srl=catalog_srl_initial.copy()#catalog_srl_initial[['datetime','source','MW']].copy()
catalog_srl=catalog_srl.drop(columns=['Date','Time','EventID'])

# Merge, sort and reset index
catalog_merged=pandas.concat([catalog_es,catalog_srl],ignore_index=True) #Concatenate the catalogs and ignore the index
catalog_merged=catalog_merged.sort_values(by='datetime') # Sort values based in date and time

catalog_merged.reset_index(drop=True,inplace=True) # Reset index after sorting
###############################################################################


#%% Part 2 ###################################################################
catalog_merged['double_flag']=False #Label for events identified by both methods

# The criterium for considering two consecutives events as the same is that they delay is =< 1.5 seconds
for index in range(1,len(catalog_merged)):

    # New condition: delay between events must be lower that 2 seconds and must be from different catalogs
    if ((catalog_merged.datetime[index] - catalog_merged.datetime[index-1]) < pandas.Timedelta('2 second')
        and (catalog_merged.source[index] != catalog_merged.source[index-1])):
        catalog_merged.at[index,'double_flag'] = True
        catalog_merged.at[index-1,'double_flag'] = True

# Sppliting the merged catalo into two new catalog. One of single, and other of double flagged
catalog_merged_single=catalog_merged[catalog_merged['double_flag']==False]
catalog_merged_single.reset_index(drop=True,inplace=True)

catalog_merged_double=catalog_merged[catalog_merged['double_flag']==True]
catalog_merged_double.reset_index(drop=True,inplace=True)


#%% Part 3 ###################################################################

### Header ###################

# This function gets a datetime string and outputs the same of 1 or 2 files containing the waveforms X-seconds around the datetime
from aux_modules import time2file_name, load_and_slice, plot_3c

# Import catalog of files names and times
catalog_time_n_files = pandas.read_csv("catalog_time_n_files.csv")
catalog_time_n_files['date_time']= list(map(lambda x: obspy.core.utcdatetime.UTCDateTime(x), catalog_time_n_files['date_time'].values))

# Folder containing the SEGY files
folder_segy_files='/Volumes/toc2me/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C_SEG2/'

# Folder to store the plots
folder_plots='plots_double/'

#%% This is a test to confirm that the double-events happens in pairs and not in trio
# The test was conducted and the hypothesis of events in pairs only was confirmed.
# This will make the following analysis easier
for index_double in range(2,len(catalog_merged_double)):
    if (((catalog_merged_double.datetime[index_double] - catalog_merged_double.datetime[index_double-1]) < pandas.Timedelta('2 second'))
        and ((catalog_merged_double.datetime[index_double] - catalog_merged_double.datetime[index_double-2]) < pandas.Timedelta('2 second'))):
        print('annomalie at',index_double)

#%% Exporting the catalogs for future use
catalog_merged_double.to_csv('catalog_merged_double.csv',sep=',',line_terminator='\n',index=True)
catalog_merged_single.to_csv('catalog_merged_single.csv',sep=',',line_terminator='\n',index=True)


#%% Input datetime and get file or files 2 load
for event_index in tqdm(range(int(len(catalog_merged_double)/2))):
#    segy_files=time2file_name(input_time=catalog_merged_double.datetime[2*event_index],sec_before=3,sec_after=3,catalog_time_n_files=catalog_time_n_files)
    #print('------------Using slots ',2*event_index,2*event_index+1)
    # get the one or 2 files containing the waveforms around the time of interest
    segy_files=time2file_name(input_time=catalog_merged_double.datetime[2*event_index],sec_before=3,sec_after=3,catalog_time_n_files=catalog_time_n_files)
    #print(segy_files,"\n")
    # Use the function to load and slice around the time of interest
    ms_data=load_and_slice(files2load=segy_files,folder=folder_segy_files,datetime2load=UTCDateTime(catalog_merged_double.datetime[2*event_index]),sec_before=3,sec_after=3)
    #print('Start and endtime',ms_data[0].stats.starttime,ms_data[0].stats.endtime,"\n")


    # Double checking which datetime is from SRL and which is from ES
    if catalog_merged_double.source[2*event_index]=='srl':
        time_srl=numpy.datetime64(catalog_merged_double.datetime[2*event_index])
    elif catalog_merged_double.source[2*event_index]=='es':
        time_es=numpy.datetime64(catalog_merged_double.datetime[2*event_index])

    if catalog_merged_double.source[2*event_index+1]=='srl':
        time_srl=numpy.datetime64(catalog_merged_double.datetime[2*event_index+1])
    elif catalog_merged_double.source[2*event_index+1]=='es':
        time_es=numpy.datetime64(catalog_merged_double.datetime[2*event_index+1])

    #print('time es and srl',time_es,time_srl,"\n")

    # Plot the waveforms, times for ES and SRL and energy stack
    plot_3c(ms_data=ms_data,time_es=time_es,time_srl=time_srl,output_folder=folder_plots)

#%%
