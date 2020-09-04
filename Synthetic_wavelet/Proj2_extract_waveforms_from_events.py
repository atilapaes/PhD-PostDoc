#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 20:13:59 2020

@author: atilapaes

Extract waveform from real events(17 sample events)
"""
import obspy, pandas, numpy
from tqdm import tqdm

folder_ms_files='Toc2me_events/'
event_data=pandas.read_csv('Proj2_event_data.csv')

# Time in sec to slice the waveform after the picking time
sec_range=0.1

# Array to store wavelet data
wavelets=numpy.zeros((17,51))


for index_event in tqdm(range(17)):

    # Read the SEGY files
    ms_data=obspy.read(folder_ms_files+event_data.file_name[index_event])


    # Slice the event from picking to sec_range(0.1) secs after the P-picking
    slice_p=ms_data[event_data.gph_first_arrival[index_event]].slice(starttime=ms_data[0].stats.starttime+event_data.atp_p_sec[index_event],
                                                                     endtime=ms_data[0].stats.starttime+event_data.atp_p_sec[index_event]+sec_range)#.data
    # Plot the sliced waveform (for visual quality control)
    slice_p.plot()

    # Store the data from each waveform to the numpy array
    wavelets[index_event][:]=slice_p.data


# Save data in external file
numpy.save('Proj2_wavelets.npy',wavelets)
