#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 18:27:50 2020

@author: atilapaes

Merging catalogs V2
This script is intended to provide a more reliable method to merge two different catalogs.
Project 1 is intended to merge two different catalogs, and then split it in two subcatalogs.
One of single events and the other of double events (events up to 2 seconds apart from its neighbour)


Part 1 - Load and treatment of thw two catalogs. Similar to the previously developed Proj01_catalog_merger_double_with_MW.py)

Part 2 - Flag doubles and split the merged catalogs into single and doubles
    
Part 3 - Plot the waveforms and time for the events, export the figures for manual inspection

Part 4 - Exporting the catalogs
"""

import pandas, obspy, numpy
from obspy.core import UTCDateTime
from aux_modules import time2file_name, load_and_slice, plot_3cV2
from tqdm import tqdm

#%%### Part 1 #################################################################
# Catalog MFA and treatment
catalog_MFA=pandas.read_csv('catalog_MFA.csv') # Read file
catalog_MFA['datetime']=catalog_MFA['Date']+'T'+catalog_MFA['Time'] # Conpose a datetime format readable by pandas
catalog_MFA['datetime']=pandas.to_datetime(catalog_MFA['datetime']) #Convert the strings generated before in the datetime format used by Pandas
catalog_MFA['source']='MFA' # Creating a column to be label that will identify the sorce identification method
catalog_MFA=catalog_MFA.drop(columns=['Time','EventID']) # Treatment of MFA's catalog. The 3D coordinates and the MW are preserved

 
# Catalog ES and treatment
catalog_es=pandas.read_csv('catalog_eslsnr.csv') # Read file
catalog_es['Date']=list(map(lambda x: x.split(' ')[0], catalog_es['datetime'].values)) #Remoce microsecond information from the ES catalog
catalog_es['datetime']=list(map(lambda x: x.split('.')[0], catalog_es['datetime'].values)) #Remoce microsecond information from the ES catalog
catalog_es=catalog_es.drop(columns=['ind','quality','method'])
catalog_es['datetime']=pandas.to_datetime(catalog_es['datetime']) # Convert column with string of datetime to the format of Pandas datetime
catalog_es['source']='ES' # Creating a column to be label that will identify the sorce identification method
catalog_es.drop_duplicates(subset='datetime',inplace=True)# Treatment of ES's catalog (there was some ducplicates)
catalog_es.sort_values(by='datetime')#.reset_index()

# Merge both catalogues, sort and reset index
catalog_merged=pandas.concat([catalog_es,catalog_MFA],ignore_index=True) #Concatenate the catalogs and ignore the index
catalog_merged=catalog_merged.sort_values(by='datetime') # Sort values based in date and time
catalog_merged.reset_index(drop=True,inplace=True) # Reset index after sorting

###############################################################################
#%% Part 2 ###################################################################
catalog_merged['double_flag']=False #Label for events identified by both methods

# The criterium for considering two consecutives events as the same is that their delay is =< 3 seconds
for index in range(1,len(catalog_merged)):
       
    if ((catalog_merged.datetime[index] - catalog_merged.datetime[index-1]) < pandas.Timedelta('3 second')
        and (catalog_merged.source[index] != catalog_merged.source[index-1])):       
        catalog_merged.at[index,'double_flag'] = True
        catalog_merged.at[index-1,'double_flag'] = True

# Spliting the merged catalog into two new catalogues. One of single, and other of double flagged
catalog_merged_single=catalog_merged[catalog_merged['double_flag']==False]
catalog_merged_single.reset_index(drop=True,inplace=True)
catalog_merged_single=catalog_merged_single.drop(columns=['double_flag'])

catalog_merged_double=catalog_merged[catalog_merged['double_flag']==True]
catalog_merged_double.reset_index(drop=True,inplace=True)


#%% Part 3 ###################################################################
  
### Header ###################
# Folder containing the SEGY files
folder_segy_files='/Volumes/toc2me/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C_SEG2/'

# Folder to store the plots
folder_plots='plots_double/'

# Import catalog with SEG2 files_names and Datetimes
catalog_time_n_files = pandas.read_csv("catalog_time_n_files.csv")
catalog_time_n_files['date_time']= list(map(lambda x: obspy.core.utcdatetime.UTCDateTime(x), catalog_time_n_files['date_time'].values))

# Flagging the set of events closed to each other (considered as doubles) by a new index and store such index in the double_block column
catalog_merged_double['double_block']=0

# Initialization of the index of double events
double_index=0

# Loop for labeling the blocks os doubles (Note: some blocks contain more than 2 entries for a close Datetime)
for index_event in range(1,len(catalog_merged_double)):
    if (catalog_merged_double.datetime[index_event] - catalog_merged_double.datetime[index_event-1]) < pandas.Timedelta('2 second'):
        catalog_merged_double.double_block[index_event] = double_index
    else:
        double_index +=1
        catalog_merged_double.double_block[index_event] = double_index

# Getting the index of the blocks with 3 or more doubles
multiple_double_list=[]
index_last_block=catalog_merged_double.double_block[len(catalog_merged_double.double_block)-1] # The index of the last block of doubles

# Loop to append the index_block
for index_block in range(1,index_last_block):
    n_doubles=len(catalog_merged_double[catalog_merged_double['double_block']==index_block])
    if n_doubles >2:
        #print(index_block,n_doubles)
        multiple_double_list.append(index_block)

#%%  Marking the events for manual inspection

# Column to flag events for manual inspection
catalog_merged_double['manual']=False


for index_event in range(len(catalog_merged_double)):
    if int(catalog_merged_double.double_block[index_event]) in multiple_double_list:
        catalog_merged_double.at[index_event,'manual'] = True
        #print(index_event, catalog_merged_double.double_block[index_event],catalog_merged_double.manual[index_event])

# Plot the events for manual inspection
for multiple_index in tqdm(range(len(multiple_double_list))):#tqdm(range(int(len(catalog_merged_double)/2))):
    # Get the index of the first event of a certain block with multiple (>2) doubles
    event_index=catalog_merged_double[catalog_merged_double['double_block']==multiple_double_list[multiple_index]].index.values[0]
    
    # Get the one or 2 files containing the waveforms around the time of interest
    segy_files=time2file_name(input_time=catalog_merged_double.datetime[event_index],sec_before=3,sec_after=3,catalog_time_n_files=catalog_time_n_files)
    #print(segy_files,"\n")
    
    # Use the function to load and slice around the time of interest
    ms_data=load_and_slice(files2load=segy_files,folder=folder_segy_files,datetime2load=UTCDateTime(catalog_merged_double.datetime[event_index]),sec_before=5,sec_after=5)
    #print('Start and endtime',ms_data[0].stats.starttime,ms_data[0].stats.endtime,"\n")

    # Plot the waveforms, times for ES and MFA and energy stack
    plot_3cV2(ms_data=ms_data,block_index=multiple_double_list[multiple_index],output_folder=folder_plots) 

#%% Part 4 - Treatement of double events with two entries

# The catalog
catalog_merged_double_two_events=catalog_merged_double[catalog_merged_double['manual']==False]
catalog_merged_double_two_events.reset_index(drop=True,inplace=True)

# At this point, I know that this catalog is composed by couples. Each couple has a 'ES' and 'MFA' entry of each event.
# The data treatment is 1) to preserve the info from MFA, and 2) rename the source for 'Both', 
# 3) drop unused (from now on) columns, 4) merge it it back to the single catalog

# 1) discard duplicates and preserve the info from MFA 
catalog_merged_double_two_events=catalog_merged_double_two_events[catalog_merged_double_two_events['source']=='MFA']
catalog_merged_double_two_events.reset_index(drop=True,inplace=True) # reseting the index

# 2) Change the source from 'MFA' to 'Both'
catalog_merged_double_two_events['source']='Both'

# 3) Drop unused columns 
catalog_merged_double_two_events=catalog_merged_double_two_events.drop(columns=['double_flag','double_block','manual'])

# 4) Merge it back to the single catalog
catalog_merged_single_remerged=pandas.concat([catalog_merged_single,catalog_merged_double_two_events],ignore_index=True) #Concatenate the catalogs and ignore the index
catalog_merged_single_remerged=catalog_merged_single_remerged.sort_values(by='datetime') # Sort values based in date and time
catalog_merged_single_remerged.reset_index(drop=True,inplace=True) # Reset index after sorting

#%% Part 5 

# The catalog for manual QC
catalog_merged_double_three_plus=catalog_merged_double[catalog_merged_double['manual']==True]
catalog_merged_double_three_plus.reset_index(drop=True,inplace=True)

# Event to be used in semi-automatic QC. Not used after manual QC
#catalog_manual_qc=pandas.DataFrame(data=multiple_double_list,index=numpy.arange(len(multiple_double_list)), columns=['double_block'])
#catalog_manual_qc['event2use']=''

#%% Part 6 
# Exporting the catalogs to use in project 2

catalog_merged_single_remerged.to_csv('result_catalog_merged_single.csv',sep=',',line_terminator='\n',index=True)

catalog_merged_double_three_plus.to_csv('result_catalog_merged_manual.csv',sep=',',line_terminator='\n',index=True)

#catalog_manual_qc.to_csv('result_manual_qc.csv',sep=',',line_terminator='\n',index=True)
#%%



