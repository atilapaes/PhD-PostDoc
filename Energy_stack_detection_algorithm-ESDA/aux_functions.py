#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 23:08:59 2020

@author: Dr. Atila Paes

Project: Energy Stack Detection Algorithm (ESDA) V5
This project is written (but not exclusive) to be used with the ToC2ME dataset

Module: auxiliar_function:
All functions nedded by the main module

"""
import obspy, numpy, gc
from tqdm import tqdm
import  matplotlib.pyplot as plt
from matplotlib import gridspec

#%%
import io_info


#%%
def datetime2string(date_time):
    """
    Returns a string for date time in a format that is proper for using as a file name
    example:
        from: numpy.datetime64('2016-11-01T00:00:00.000000')
        filename become: 2016-11-01T00/00/00.000000Z
        
        to: 2016-11-01T00-00-00    
    """
    return(str(date_time).split('.')[0].replace(':','-'))

#%% Catalog functions
def day2files(catalog,date_str):
    """
    This function returns a list of SEG2 files for a certain day
    
    Parameters:
    -----------
    catalog: A pandas DataFrame with the columns 
            'date' (string) 
            'file_name' (string)
    
    date_str: A string of the day to get the files. Format 'YYYY-MM-DD' 
    
    Returns:
    -------    
    list_files: an 1-D numpy ndarray of files. Each entry is a string
    """
    
    list_files=catalog[catalog['date']==date_str].file_name.values
    return(list_files)

#list_files=day2files(df,'2016-12-01')
#print(list_files[0],list_files[-1])    

#%%
def load_list_seg2(list_files,folder):
    """
    Load a list of SEGY2 files and concatenate into a single stream
    """
    #loading first file
    ms_data=obspy.read(folder+list_files[0],format='SEG2')
    
    # Load and concatenated the rest of files
    for index_file in tqdm(range(1,len(list_files))):
        ms_data_aux=obspy.read(folder+list_files[index_file])
        
        for channel in range(len(ms_data)):
            ms_data[channel] += ms_data_aux[channel]  
    return(ms_data)
#%%
def mavg(signal,samples):
    """
    This function calculates the moving average of a provided 1-C signal (array).
    """
    if numpy.sum(signal) != 0:
        signal_ma=numpy.zeros(len(signal))
        #Regular signal
        for index in range(samples//2, len(signal)-samples//2):
            signal_ma[index] = (numpy.sum(signal[(index-samples//2):(index+samples//2)]))/samples
             
        #Borders of the signal
        signal_ma[0:samples//2] = 0
        signal_ma[len(signal)-samples//2:len(signal)] = 0
        signal_ma=signal_ma/signal_ma.max()
        return (signal_ma)
    else:
        return (signal)


#%%
def cf_es(ms_data,lns, normalize,mavg_apply,mavg_samples):
    """
    Module for energy stack based on a list of selected geophones.
    Builded on top of a deep copy of first channel.
    
    lns - low noise stations
    
    if samples_mavg is 0, it does 
    """
    if normalize==True:
        ms_data=ms_data.normalize()
    
    es=ms_data[0].copy()
    es.data *= 0 # trace zeroed for storing Energy Stack info
            
    # Energy stack calculation for valid gph only
    for index_station in range(len(lns)):
        #print('======',index_station,lns[index_station],lns[index_station]+68,lns[index_station]+2*68)
        es.data = numpy.add(es.data,numpy.add(numpy.square(ms_data[lns[index_station]].data),
                                    numpy.add(numpy.square(ms_data[lns[index_station]+68].data),
                                              numpy.square(ms_data[lns[index_station]+2*68].data))))
    
    if mavg_apply==True:
        es.data=mavg(signal=es.data,samples=mavg_samples)
    es=es.normalize()
    return(es)

#%%
def cf_es_all_stations(ms_data,lns, normalize,mavg_apply,mavg_samples):
    """
    Module for energy stack based on a list of selected geophones.
    Builded on top of a deep copy of first channel.
    
    lns - low noise stations
    
    if samples_mavg is 0, it does 
    """
    if normalize==True:
        ms_data=ms_data.normalize()
    
    es=ms_data[0].copy()
    es.data *= 0 # trace zeroed for storing Energy Stack info
            
    # Energy stack calculation for valid gph only
    for index_station in range(69):
        es.data = numpy.add(es.data,numpy.add(numpy.square(ms_data[index_station     ].data),
                                    numpy.add(numpy.square(ms_data[index_station  +68].data),
                                              numpy.square(ms_data[index_station+2*68].data))))
    
    if mavg_apply==True:
        es.data=mavg(signal=es.data,samples=mavg_samples)
    es=es.normalize()
    return(es)

#%%
def stations_energy(ms_data,normalize=True,mavg_apply=True,mavg_samples=50):
    """
    Builded on top of a deep copy of z channels.
    """
    
    if normalize==True:
        ms_data=ms_data.normalize()
    
    energy=ms_data[0:69].copy()
    for index_station in range(69):
        energy[index_station].data =numpy.add(numpy.square(ms_data[index_station     ].data),
                                    numpy.add(numpy.square(ms_data[index_station  +68].data),
                                              numpy.square(ms_data[index_station+2*68].data)))
        if mavg_apply == True:
            energy[index_station].data = mavg(signal=energy[index_station].data,samples=mavg_samples)
                
    energy=energy.normalize()
    
    return(energy)

#%%
def automatic_qc(es_slice,min_peak_height,min_peak_dist,max_n_peaks,peak_list_es):
    """
    Parameters
    ----------
    es_slice :          Obspy trace about 8 seconds long
    min_peak_height:    Decimal threshold for triggering a peak.
    min_peak_dist:      Minimum distance (in samples) for triggering a second peak 
    max_n_peaks:        Maximum number of peaks to identify the ES characteristic double-hump
        
    Returns
    -------
    test_qc:   logic
                If confirmed the event is a true positive -> returns True, otherwise returns False

    """
    
    # Normalize before looking for peak (local normalization)
    es_slice=es_slice.normalize()
    peak_list_es=detect_peaks(x=es_slice.data,mph=min_peak_height,mpd=min_peak_dist, show=False)
                
    
    #### QC 1 ####################
    # Test  the number of peaks in ES must be 3 or less
    if len(peak_list_es) <= max_n_peaks:
        test_qc1=True
    else:
        test_qc1=False
        
    #### QC 2 ####################
    # Test Considering a sample of 8 seconds, there are 4k samples.
    # and it is expected that the peaks are between 2 and 6 seconds. 
    
    # Count number of peaks where this conditon is valid
    count_peak_into_range=0
    for index_peak in range(len(peak_list_es)):
        if (peak_list_es[index_peak] >=1000) and (peak_list_es[index_peak] <=3000):
            count_peak_into_range +=1
    
    # The test in QC 2 is valid if at least 2  peaks are into this interval
    if count_peak_into_range >=2:
        test_qc2=True
    else:
        test_qc2=False
    
    ### Decision of True Postive
    # To be considered a True positive, the potential event must satisfy both QCs
    test_qc = test_qc1 and test_qc2    

    return(test_qc)


#%%
def plot_energy_station(energy_station,es,output_folder,plot_title):
    energy_station=energy_station.normalize()
    time=numpy.arange(numpy.datetime64(energy_station[0].stats.starttime),numpy.datetime64(energy_station[0].stats.endtime)+numpy.timedelta64(int(energy_station[0].stats.delta*1000),'ms'), numpy.timedelta64(int(energy_station[0].stats.delta*1000),'ms'))    
    
    fig = plt.figure(1,figsize=(6,8))
    fig.subplots_adjust(left=0.08,right=0.95, bottom=0.05,top=0.95,wspace=0.1, hspace=0.10)
    gs = gridspec.GridSpec(2, 1, height_ratios=[7, 1]) 
    
    ax1=plt.subplot(gs[0])#(2,1,1)
    for gph_index in range(69):               
        
        if gph_index in io_info.lns:
            ax1.plot_date(time,gph_index+0.45*energy_station[gph_index].data,'k',lw=0.3)
        else:
            ax1.plot_date(time,gph_index+0.45*energy_station[gph_index].data,'r',lw=0.3)
            
    ax1.set_yticks(numpy.arange(0,71,2))
    ax1.tick_params(axis='both',labelsize='xx-small')
    ax1.set_title(plot_title,fontsize='small')    
    ax1.set_ylabel('Stations energy',fontsize='small')
    ax1.set_ylim(-1,70)
    ax1.set_xlim(time[0],time[-1])
    
    ax2=plt.subplot(gs[1])
    
    # Calculating the energy stack of the sample
    ax2.plot(time,es.data,'k',lw=1,label='Energy stack')
    ax2.legend(loc='upper left')
    #ax2.set_xticklabels(fontsize='x-small')
    ax2.set_xlim(time[0],time[-1])
    ax2.tick_params(axis='both', which='major', labelsize=6)
    
    plt.savefig(output_folder +plot_title+'.jpg',dpi=300,format='jpg')
    plt.close()
    
    # Clear the memory after saving figure
    fig.clf()
    plt.close('all')
    gc.collect()
    
    #print('Plot done')
    
    return()


#%%

def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    import numpy

    """Detect peaks in data based on their amplitude and other features.


    Detect peaks in data based on their amplitude and other features.
    http://nbviewer.jupyter.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Atila's example
    u=detect_peaks.detect_peaks(data1, mph=numpy.mean(data1)+1.5*numpy.std(data1), mpd=100, show=True)


    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`

    The function can handle NaN's

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    """

    x = numpy.atleast_1d(x).astype('float64')
    if x.size < 3:
        return numpy.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = numpy.where(numpy.isnan(x))[0]
    if indnan.size:
        x[indnan] = numpy.inf
        dx[numpy.where(numpy.isnan(dx))[0]] = numpy.inf
    ine, ire, ife = numpy.array([[], [], []], dtype=int)
    if not edge:
        ine = numpy.where((numpy.hstack((dx, 0)) < 0) & (numpy.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = numpy.where((numpy.hstack((dx, 0)) <= 0) & (numpy.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = numpy.where((numpy.hstack((dx, 0)) < 0) & (numpy.hstack((0, dx)) >= 0))[0]
    ind = numpy.unique(numpy.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[numpy.in1d(ind, numpy.unique(numpy.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = numpy.min(numpy.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = numpy.delete(ind, numpy.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[numpy.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = numpy.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = numpy.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = numpy.nan
        if valley:
            x = -x
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind


def _plot(x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    import numpy
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(x, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(ind, x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.axhline(mph,color='cyan')

            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ax.set_xlim(-.02*x.size, x.size*1.02-1)
        ymin, ymax = x[numpy.isfinite(x)].min(), x[numpy.isfinite(x)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1*yrange, ymax + 0.1*yrange)
        ax.set_xlabel('data #', fontsize=14)
        ax.set_ylabel('Amplitude', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()
#%%