#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 13:11:15 2020

@author: atilapaes

This scrip reads the waveforms from previously generated Proj2_wavelets, 
calculates the best fit with a Berlage waveform. Plots the event spectrum and the fitting

"""
import numpy
import matplotlib.pyplot as plt
from Berlage_wavelet import berlage2
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab
#%%
index_event=0

# Import the .npy file generated before
wavelet=numpy.load('Proj2_wavelets.npy')

#Array for time in seconds
t=numpy.arange(0,0.102,0.002)

#%% Getting the frequency of maximum amplitude

# Data taper
wl_tap=numpy.concatenate((numpy.zeros(500),wavelet[index_event]),axis=0)

spec, freqs = mlab.magnitude_spectrum(x=wl_tap, Fs=1/0.002)
spec /=numpy.max(spec)
 
freq_max_amp=freqs[numpy.argmax(spec)]

#%% setting the initial values
init_vals= [freq_max_amp,float(-100000),float(2),float(50),-numpy.pi/2]

#%% Calculating the parameter for the best fit
best_vals, covar = curve_fit(berlage2,t,wavelet[index_event],p0=init_vals)

#%% Calculating the Berlage wavelet with the parameters for the best fit
berlage2_best=berlage2(t,best_vals[0],best_vals[1],best_vals[2],best_vals[3],best_vals[4])

#%% The plot
plt.figure(3,figsize=(6,6))
plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95,wspace=0.3,hspace=0.3)

ax1=plt.subplot(2,1,1)
ax1.set_title('Event '+str(index_event))
ax1.plot(freqs,spec,'-o',lw=0.5,ms=2)
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Magnitude (normalized)')

ax2=plt.subplot(2,1,2)
#ax2.set_title('Event '+str(index_event))
ax2.plot(t,wavelet[index_event],'-ob',lw=0.5,ms=2,label='Original data')
ax2.plot(t,berlage2_best,'-or',lw=0.5,ms=2,label='Best fit')    
ax2.legend(loc='upper left')
ax2.set_xlim(-0.005,0.15)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Amplitude (cte *nm/s)')

col_labels=['Initial','Best']
row_labels=['Freq','Ampl.','n','alpha','phi']
table_vals=[[str("%.2f" % init_vals[0]),str("%.2f" % best_vals[0])],
            [str("%.2f" % init_vals[1]),str("%.2f" % best_vals[1])],
            [str("%.2f" % init_vals[2]),str("%.2f" % best_vals[2])],
            [str("%.2f" % init_vals[3]),str("%.2f" % best_vals[3])],
            [str("%.2f" % init_vals[4]),str("%.2f" % best_vals[4])]]

the_table = plt.table(cellText=table_vals,
                  colWidths = [0.1]*3,
                  rowLabels=row_labels,
                  colLabels=col_labels,
                  loc='upper right').scale(1.2, 1.2) 
plt.savefig('event'+str(index_event)+'.png',format='png',dpi=300)
plt.show()
plt.close()
