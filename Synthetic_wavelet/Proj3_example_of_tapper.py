#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 21:07:52 2020

@author: atilapaes
"""

import numpy
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#%%
index_event=4

# Import the .npy file generated before
wavelet=numpy.load('Proj2_wavelets.npy')

#Array for time in seconds
t=numpy.arange(0,0.102,0.002)

#%% Original data
wl=wavelet[index_event]
spec, freqs = mlab.magnitude_spectrum(x=wl, Fs=1/0.002)
spec /=numpy.max(spec) # normalizing for comparative plot
freq_max_amp=freqs[numpy.argmax(spec)]

#%% Tapped data
wl_tap=numpy.concatenate((numpy.zeros(500),wavelet[index_event]),axis=0)
spec_tap, freqs_tap = mlab.magnitude_spectrum(x=wl_tap, Fs=1/0.002)
spec_tap /=numpy.max(spec_tap) #normalizing for comparative plot
freq_max_amp_tap=freqs_tap[numpy.argmax(spec_tap)]
#%%

plt.figure(1)
plt.title('Event '+str(index_event))
plt.plot(freqs_tap,spec_tap,'-or',lw=0.5,ms=2,label='Taper data')
plt.plot(freqs,spec,'-sb',lw=0.5,ms=4,label='Original data')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Normalized mag spectrum')
plt.legend(loc='upper right')
plt.show()
plt.close()