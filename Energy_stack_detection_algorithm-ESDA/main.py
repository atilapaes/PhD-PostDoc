#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 23:06:49 2020

@author: Dr. Atila Paes

Project: Energy Stack Detection Algorithm (ESDA) V5
This project is written (but not exclusive) to be used with the ToC2ME dataset


Module: main:
This is the main module to use  Energy Stack Detection Algorithm (ESDA) 
using the workflows developed.

"""

import pandas, obspy, numpy
import matplotlib.pyplot as plt
#%%
import io_info, aux_functions

#%%
def run_detection(list_days2process):#(day2process):
    
    # Loop over the list of days
    for day2process in list_days2process:
        
        print('==================================')
        print('=== Processing day',day2process)
        print('==================================')
        
        #Read catalog of files
        catalog_seg2_files=pandas.read_csv(io_info.catalog_seg2_files)
        
        #Get the list for a certain file
        list_files=aux_functions.day2files(catalog_seg2_files,day2process)
        
        # Divide the files of a day in chucks to save memory
        list_files_slice=numpy.array_split(list_files,io_info.slice_day)
        
        #%% Data loading and treatment
        
        # Loop over the number of slices a day is divided in
        for index_slice in range(len(list_files_slice)):
            print('==================================')
            print('=== Processing slice',index_slice,' of ',len(list_files_slice))
            print('==================================')
            
            print('== Loading SEG2 files...')
            ms_data=aux_functions.load_list_seg2(list_files=list_files_slice[index_slice],folder=io_info.folder_segy2_files)
            print('== Loading done')
            
            #%%
            print('==== Starting data treatment')
            print('== Filtering')
            ms_data=ms_data.filter('bandpass', freqmin=10, freqmax=80, corners=4, zerophase=True)
            
            print('== Stack of energy') # Using mavg on es cf makes the curve less crisp and decreases the number of false positives
            #es=aux_functions.cf_es(ms_data,lns=io_info.lns, normalize=True,mavg_apply=True,mavg_samples=io_info.mavg_samples)
            
            es_all_stations=aux_functions.cf_es_all_stations(ms_data,lns=io_info.lns, normalize=True,mavg_apply=True,mavg_samples=io_info.mavg_samples)
            #%%
            
            print('== Stations energy')
            energy_stations=aux_functions.stations_energy(ms_data=ms_data,normalize=True,mavg_apply=False,mavg_samples=io_info.mavg_samples)
            
            
            # TIME - Using numpy datetime64 - When needed for obspy, convert using obspy lib
            time=numpy.arange(numpy.datetime64(ms_data[0].stats.starttime),numpy.datetime64(ms_data[0].stats.endtime)+numpy.timedelta64(2,'ms'), numpy.timedelta64(2,'ms'))
                
            #print('== Data treatment ')
            #energy.plot(inplace=False)
            #plt.plot(es,lw=0.1)
            
            #%%
            # ES analysis - Threshold chose as % of maximum for 1 hour, interactive for each hour
            print('=== Starting detection')
            print('== Peak detection')
            mph=io_info.mph
            mpd=int(io_info.peak_distance)
                
            peak_list=aux_functions.detect_peaks(x=es_all_stations.data,mph=mph,mpd=mpd, show=False)
            
            #%% With the peak list and the time vector, I ca slice the energy and 
            # 1) save the plots
            # 2) Automatic QC
            
            print('== Automatic QC')
            events_detected=pandas.DataFrame(data=None, columns=['datetime','true_positive'])
            
            for index_pe in range(len(peak_list)):
                
                # Slice the energy per station and the stack os energy
                energy_stations_slice=energy_stations.slice(starttime=obspy.core.utcdatetime.UTCDateTime(str(time[peak_list[index_pe]])) - io_info.slice_sec_b,
                                                                                  endtime=obspy.core.utcdatetime.UTCDateTime(str(time[peak_list[index_pe]])) + io_info.slice_sec_a)         
                
                es_slice=es_all_stations.slice(starttime=obspy.core.utcdatetime.UTCDateTime(str(time[peak_list[index_pe]])) - io_info.slice_sec_b,
                                                                                  endtime=obspy.core.utcdatetime.UTCDateTime(str(time[peak_list[index_pe]])) + io_info.slice_sec_a)         
                                    
                # Normalize before looking for peak (local normalization)
                es_slice=es_slice.normalize()
                peak_list_es=aux_functions.detect_peaks(x=es_slice.data,mph=0.3,mpd=300, show=False)
                
                # Function for automatic qualitu control
                true_positive_event=aux_functions.automatic_qc(es_slice=es_slice,min_peak_height=io_info.min_peak_height,min_peak_dist=io_info.min_peak_dist,max_n_peaks=io_info.max_n_peaks,peak_list_es=peak_list_es)
                
                # If True positive, plot event and save its infor into the catalog
                if true_positive_event == True:
                    
                    # Plot the energy stack
                    aux_functions.plot_energy_station(energy_station=energy_stations_slice,es=es_slice,output_folder=io_info.folder_pe,plot_title=aux_functions.datetime2string(time[peak_list[index_pe]]))
            
                    # Append the event
                    events_detected=events_detected.append({'datetime':time[peak_list[index_pe]],'true_positive':''},ignore_index=True)
            
            print('== Exporting results')
            # Save the dataframe with Potential events for a slice of a day
            events_detected.to_csv(io_info.folder_catalog_events + day2process + '_slice_'+str(index_slice) + '_manual_qc.csv',sep=',',line_terminator='\n',index=False)    
    
    print('==== Processing done ====')        
#%%
run_detection(list_days2process=['2016-11-02',
                                '2016-11-03',
                                '2016-11-04',
                                '2016-11-05',
                                '2016-11-06',
                                '2016-11-07',
                                '2016-11-08',
                                '2016-11-09',
                                '2016-11-10',
                                '2016-11-11',
                                '2016-11-12',
                                '2016-11-13',
                                '2016-11-14',
                                '2016-11-15',
                                '2016-11-16',
                                '2016-11-17',
                                '2016-11-18',
                                '2016-11-19',
                                '2016-11-20',
                                '2016-11-21',
                                '2016-11-22',
                                '2016-11-23',
                                '2016-11-24',
                                '2016-11-25',
                                '2016-11-26',
                                '2016-11-27',
                                '2016-11-28',
                                '2016-11-29',
                                '2016-11-30'])