#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 21:09:47 2019

@author: atilapaes

Calculate the SNR of the each Z channel
builded primarly for the 3C gph, channel Z

This version loads the ms_data, calculates the sum of squares elements and 
store in the memory only the sum of squared values and number of valid elements
"""

import numpy, obspy, pandas
from tqdm import tqdm
import matplotlib.pyplot as plt

#%% Parameters
folder_3c='/Volumes/toc2me/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C_SEG2/'
#day2process=21

#%% SPECIFIC FUNTIONS
def trace_rms(data):
    """
    Data: 1C numpy array that can contanin "None" values
    """    
    rms=numpy.sqrt(numpy.nansum(numpy.square(data))/len(data))
    return(rms)
    
###############################################################################    
### Load and concatenate list of files   ###################################### 
def load_list_files(folder_3c, list_3c):
    """
    Load and concatenate list of files in a one day-long catalog
    Return just Z component    
    
    """
    print('Loading ',folder_3c+list_3c['file_name'][0])
    ms_data=obspy.read(folder_3c+list_3c['file_name'][0])
    ms_data=ms_data[0:69]
    
    for file_index in tqdm(range(1,len(list_3c))):
        #print('')
        #print('Loading ',list_3c['file_name'][file_index])
        ms_data_aux=obspy.read(folder_3c+list_3c['file_name'][file_index])
        ms_data_aux=ms_data_aux[0:69]
        
        for ch_index in range(69):
            ms_data[ch_index] += ms_data_aux[ch_index]
            
    return(ms_data)    
###############################################################################

#%% Catalogs readind and treatment

#### Eaton's catalog of events and its treatment
c_eaton=pandas.read_csv('Eaton2018.csv')
c_eaton['datetime']=c_eaton['Date']+'T'+c_eaton['Time']
c_eaton['datetime']=pandas.to_datetime(c_eaton['datetime'])


#### Atila's catalog of events and its treatment
c_atila1=pandas.read_csv('catalog_atila1.csv')
#Make a column of date in str format
c_atila1['date_str']= list(map(lambda x: x.split(' ')[0], c_atila1['datetime'].values))

c_atila1['datetime']=pandas.to_datetime(c_atila1['datetime'])
c_atila1['date']=c_atila1['datetime'].dt.date 
############################################################################

#### 3C file name of day list

### Catalog with catalogs of file names
catalog_days=pandas.read_csv('file_list/catalog_list_v2.csv')

### Loop over a day
for day_index in range(15):#len(catalog_days)
    
    ### Load list of 1 minute files to load
    list_3c=pandas.read_csv('file_list/'+catalog_days['gph3c'][day_index])

    ### List of events to exclude data from the waveforms
    #c_eaton_day=c_eaton[c_eaton['Date']==catalog_days['day_label'][day_index]]
    #c_atila1_day=c_atila1[c_atila1['date_str']==catalog_days['day_label'][day_index]]
    
    ## One rms per day (all files included) 
    sum_square=numpy.zeros(69)
    n_files_used=numpy.zeros(69)

    #### LOAD FILES AND PARTIAL CALCULATION ###################################
    for file_index in tqdm(range(len(list_3c))):#
        ms_data=obspy.read(folder_3c+list_3c['file_name'][file_index])
        ms_data=ms_data[0:69]
        
        #%% Do not use data from files with event inside
        # Start and end time of the ms data file
        s_time=pandas.Timestamp(str(ms_data[0].stats.starttime).split('.')[0])
        e_time=pandas.Timestamp(str(ms_data[0].stats.endtime).split('.')[0])

        # Slice both catalogs to check if there is at least 1 event into the file
        c_eaton_sliced=c_eaton[(c_eaton['datetime'] >= s_time) & (c_eaton['datetime'] < e_time)]
        c_atila1_sliced=c_atila1[(c_atila1['datetime']>=s_time) & (c_atila1['datetime'] < e_time)]

        if (len(c_eaton_sliced) == 0) and (len(c_atila1_sliced) == 0):
        
            ### Calculate and store the sum of squares amplitudes and number o elements
        
            for ch_index in range(69):
                # Verify if the trace is not null
                if numpy.nansum(ms_data[ch_index].data) != 0:
                    sum_square[ch_index] += numpy.nansum(numpy.square(ms_data[ch_index].data))
                    n_files_used[ch_index] += 1
    
    ##########################################################################
    
    
    ### Calculate RMS after accumulating info from a whole day ################
    rms=numpy.zeros(69)
    for ch_index in range(69):
        if n_files_used[ch_index] !=0:
            rms[ch_index] =  numpy.sqrt(sum_square[ch_index]/(30000*n_files_used[ch_index]))
        else:
            rms[ch_index]=0
            
    ### Merging Numpy arrays in a DF to export
    day_snr_df=pandas.DataFrame({'rms':rms,'sum_square':sum_square, 'n_files_used':n_files_used},index=numpy.arange(0,69,1))
            
    ### Export RMS for the specific day    
    day_snr_df.to_csv('results/' + catalog_days['day_label'][day_index] + '-3C-RMS-raw.csv',sep=',',line_terminator='\n',index=True)        

#%%



    