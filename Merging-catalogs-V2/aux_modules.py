#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 21:02:14 2020

@author: atilapaes
"""
import obspy, numpy, gc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec

#%%
def time2file_name(input_time,sec_before,sec_after,catalog_time_n_files):
    """
    This function get a input_time and returns which file(s) contain(s) the 
    time window from [input_time - sec_before, input_time - sec_after]
    
    INPUTS:
    input_time: pandas.tslib.Timestamp or Obspy datetime
    sec_before: number of seconds before input time (int)
    sec_after: number of seconds after input time (int)
    catalog_time_n_files: contains 2 columns: file_name(str), date_time(MUST BE pandas datetime type)
    
    OUTPUT:
    files2load: list of 1 or 2 file names to load
    """
    
    # Step 1: Get catalog index for input_time.
    catalog_index=catalog_time_n_files[(catalog_time_n_files['date_time'] <= input_time) &
                                       (catalog_time_n_files['date_time'] + 59.998 >= input_time) ].index[0]
    
    # STEP 2: Get file name of specified index and a secondary file (before or after) if necessary
    files2load=[catalog_time_n_files.file_name[catalog_index]]
    if input_time.second-sec_before < 0:
        files2load.append(catalog_time_n_files.file_name[catalog_index-1])
    elif input_time.second+sec_after > 60:
        files2load.append(catalog_time_n_files.file_name[catalog_index+1])
    
    return(files2load)

#############################################################################

def load_and_slice(files2load,folder,datetime2load,sec_before,sec_after):
    """
    This function load 1 or 2 files, concatenate them in case of 2 files, and slice around the datetime of interest
    """    
    # Loading first file
    ms_data=obspy.read(folder+files2load[0])

    # In case there are two files, load and concatenate them
    if len(files2load)==2:
        ms_data_aux=obspy.read(folder+files2load[1])
        for index_ch in range(len(ms_data)):
            ms_data[index_ch] += ms_data_aux[index_ch]
    
    # Slice around time of interest
    ms_data=ms_data.slice(starttime=datetime2load - sec_before, endtime=datetime2load + sec_after)
    
    return(ms_data)

#%% #############################################################################
def cf_es(ms_data):
    """
    Returns the Energy stack trace using the normalized stream
    
    """
    
    es=ms_data[0].copy()
    es.data *=0
    ms_data = ms_data.normalize()    
    for index in range(len(ms_data)):
        es.data += numpy.square(ms_data[index].data)
    return(es)

#%% #############################################################################
def plot_3c(ms_data,time_es,time_srl,output_folder):
    ms_data=ms_data.normalize()
    print('len ms_data',len(ms_data[0].data))
    time=numpy.arange(numpy.datetime64(ms_data[0].stats.starttime),numpy.datetime64(ms_data[0].stats.endtime)+numpy.timedelta64(int(ms_data[0].stats.delta*1000),'ms'), numpy.timedelta64(int(ms_data[0].stats.delta*1000),'ms'))    
    print('len time',len(time))
    
    fig = plt.figure(1,figsize=(8,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 1]) 
    
    ax1=plt.subplot(gs[0])#(2,1,1)
    for gph_index in range(69):               
        
        # This if is used to use only a single legend of each channel
        if gph_index == 0:
            ax1.plot(time,gph_index+0.45*ms_data[gph_index].data,'r',lw=0.1,label='Z')
            ax1.plot(time,gph_index+0.45*ms_data[gph_index+69].data,'g',lw=0.1,label='H1')
            ax1.plot(time,gph_index+0.45*ms_data[gph_index+2*69].data,'b',lw=0.1,label='H2')
        else:
            ax1.plot(time,gph_index+0.45*ms_data[gph_index].data,'r',lw=0.1)
            ax1.plot(time,gph_index+0.45*ms_data[gph_index+69].data,'g',lw=0.1)
            ax1.plot(time,gph_index+0.45*ms_data[gph_index+2*69].data,'b',lw=0.1)
        
    ax1.set_yticks(numpy.arange(0,71,2))
    ax1.tick_params(axis='both',labelsize='xx-small')
    #ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m:%S'))
    ax1.set_title(str(time_srl))    
    
    
    ax1.axvline(x=time_es,c='r',lw=2,label='ES')
    ax1.axvline(x=time_srl,c='b',lw=2,label='srl')
    
    ax1.legend(loc='upper left')
    ax1.yaxis.tick_right()    
    ax1.set_ylim(-1,70)
    
    ax2=plt.subplot(gs[1])
    
    # Calculating the energy stack of the sample
    es=cf_es(ms_data)
    ax2.plot(time,es.data,'k',lw=1,label='Energy stack')
    ax2.legend(loc='upper left')
    
    plt.savefig(output_folder +str(time_srl)+'.jpg',dpi=100,format='jpg')
    plt.close()
    
    # Clear the memory after saving figure
    fig.clf()
    plt.close('all')
    gc.collect()
    
    return()

#plot_3c(ms_data=ms_data,time_es=time_es,time_srl=time_srl,output_folder=folder_plots)        
#############################################################################