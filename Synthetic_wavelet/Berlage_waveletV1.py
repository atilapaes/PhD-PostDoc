#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:44:25 2020

@author: atilapaes

"""
import numpy
import matplotlib.pyplot as plt
#from matplotlib.patches import Rectangle
#import matplotlib.pylab as pylab


def berlage(f,length,dt,A,n,alpha,phi,plot=True):
    """
    This function generates A Berlage wavelet
    It returns two vectors. The time and the Berlage wavelet


    It models the function: 
        A*H(t)*T^n*exp(-alpha*t)*cos(2.Pi*f0*t+phi)
    
    -----------------------------------------------
    f   : float
        Oscillation frequency
    
    length  : int
        vector length in seconds
        Will be centered aroung zero
        
    dt      : float
        time sampling
    
    A       : float
        wave amplitude
    
    n       : float (but usually a positive constant)
            
    
    alpha   : float (nonnegative real cte)
            Decay factor
    
    phi     : float 
            phase

    Returns
    -------
    vector : (N,) ndarray
        Array of length `points` in shape of ricker curve.
    """
    
    t=numpy.arange(-length/2,length/2,dt)
    H=numpy.heaviside(t,0)
    w=numpy.zeros(int(length/dt))
    
    for index_t in range(len(w)):
        w[index_t]=A*H[index_t]*((t[index_t])**n)*numpy.exp(-alpha*t[index_t])*numpy.cos(2*numpy.pi*f*t[index_t]+phi)
    
    
    if plot==True:
        plt.plot(t,w,'.-')
    return t, w   
    
t,w=berlage(f=50,length=0.1,dt=0.002,A=100000,n=2,alpha=100,phi=-numpy.pi/2)   
    