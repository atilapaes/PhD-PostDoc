#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 22:25:10 2020

@author: Dr. Atila Paes

Project: Energy Stack Detection Algorithm (ESDA) V5
This project is written (but not exclusive) to be used with the ToC2ME dataset


Module: io_info:
This module has all information about folders and files to import for the 
Energy Stack Detection Algorithm (ESDA) using the workflows developed 
by Dr. Atila Paes.

NOTE: For using in a different computer and/or folder structure, this is the only file
that needs to be updated.

"""

#%% Processing parameters 
# How many chunks of files to divide in one day (each day has ~36 GB). 
slice_day=24 # Using 4 it takes about 6 min to load the files from USB3 and uses about 9GB RAM

#%% Folder for IO
folder_segy2_files='/Volumes/toc2me/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C/Fox_Creek_3D_2016_3C_SEG2/'
folder_output_detection='output_detection/'
folder_plots='plots/'


#%% Catalogs
catalog_seg2_files='catalogs/catalog_seg2.csv'

#%% geophones
# Low noise stations
lns=[0,1,2,3,4,5,6,8,11,14,15,16,19,21,23,24,27,29,30,34,35,37,40,41,42,45,48,49,50,51,56,57,58,59,62,66,67]

#%% processing parameters

mavg_samples=50

#%% picking parameters
peak_distance=750
mph=0.1 #input_float_number('Trigger level (decimal)? ')


slice_sec_b=4 # slice second before
slice_sec_a=4 # slice second after

folder_pe='potential_events_plots/'
folder_catalog_events='catalog_events/'

#%% QC parameters
min_peak_height=0.3
min_peak_dist=300
max_n_peaks=3