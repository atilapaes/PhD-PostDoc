#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 01:48:07 2020

@author: atilapaes
"""
import pandas

#%% Data importing and preparation

# Import catalog of single events (no duplicate present)
catalog = pandas.read_csv('result_catalog_merged_ES_MFA_FINAL.csv',index_col='Unnamed: 0')
catalog['datetime']=pandas.to_datetime(catalog['datetime'])

#%% Generate the data

# Printed using set(catalog['Date'].values)

list_days=['2016-10-26', '2016-10-27', '2016-10-28', '2016-10-29', '2016-10-30', '2016-10-31', 
 '2016-11-01', '2016-11-02', '2016-11-03','2016-11-04', '2016-11-05', '2016-11-06',
 '2016-11-07', '2016-11-08', '2016-11-09', '2016-11-10', '2016-11-11', '2016-11-12', '2016-11-13',
 '2016-11-14', '2016-11-15', '2016-11-16', '2016-11-17', '2016-11-18', '2016-11-19',
 '2016-11-20', '2016-11-21', '2016-11-22', '2016-11-23', '2016-11-24', '2016-11-25',
 '2016-11-26', '2016-11-27', '2016-11-28', '2016-11-29', '2016-11-30']

#%% Creating the dataframe
benchmark=pandas.DataFrame(data=list_days,columns=['Date'])

benchmark['ES']=''
benchmark['MFA']=''
benchmark['Both']=''

#%%
for index_day in range(len(list_days)):
    benchmark.at[index_day,'ES']=len(catalog.loc[(catalog['Date']==list_days[index_day]) & (catalog['source']=='ES')])
    benchmark.at[index_day,'MFA']=len(catalog.loc[(catalog['Date']==list_days[index_day]) & (catalog['source']=='MFA')])
    benchmark.at[index_day,'Both']=len(catalog.loc[(catalog['Date']==list_days[index_day]) & (catalog['source']=='Both')])   
    
#%% Plot
benchmark.set_index('Date').plot.bar(title='Events detected',figsize=(15,10),fontsize=12)