#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:16:31 2020

@author: atilapaes

Energy stack characteristic function development
"""
# Importing the libraries to be used
import obspy, numpy

# Reading the SEG2 file
ms_data=obspy.read('MS-FILE.dat')

# Normalize all channels to avoid signal contamination by seismic surge (see Thesis Chap 3)
ms_data=ms_data.normalize()

# Deep copy of the Z channels to use the file formart to store station energy
ms_data_squared=ms_data[0:69].copy()


#=== STEP 1 ================

# Step 1: Calculate the station energy by stacking the energy of the 3 components (Z, H1 and H2)
for index_g in range(69):
    ms_data_squared[index_g].data= numpy.square(ms_data[index_g].data) + numpy.square(ms_data[index_g+69].data) +  numpy.square(ms_data[index_g+2*69].data)


#=== STEP 2 ================

# Deep copy of a single station to use the array format to store energy stack
ms_data_es=ms_data_squared[0].copy()

# Step 2: Stack the energy of all stations
for index_g in range(1,69):
    ms_data_es.data=numpy.add(ms_data_es.data,ms_data_squared[index_g].data)


#=== STEP 3 ================

def cf_moving_avg(signal,samples=50,normalize=True):
    """
    This function calculates the moving average of a provided 1-C signal (array).
    """
    signal_ma = numpy.zeros((len(signal),))

    #Regular signal
    for index in range(samples//2, len(signal)-samples//2):
        signal_ma[index] = (numpy.sum(signal[index-samples//2:(index+samples//2)]))/samples

    #Borders of the signal are fulfilled as zero

    signal_ma[0:samples//2] = 0
    signal_ma[len(signal)-samples//2:len(signal)] = 0
    if normalize==True:
        signal_ma=signal_ma/signal_ma.max()
    return (signal_ma)

# Deep copy of a single channel of the data to use its array format
ms_data_es_mavg=ms_data_es.copy()

#Calculating the curve smoothing
ms_data_es_mavg.data=cf_moving_avg(ms_data_es_mavg.data,samples=50,normalize=True)
