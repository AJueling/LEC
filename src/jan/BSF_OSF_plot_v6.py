# Oder so starten:
#/opt/local/bin/python3.4 spectrum.py

############################################
# LOAD PACKAGES
############################################
#import os
#import scipy.io
#import scipy.stats as scst
from scipy import stats
#from functools import partial
#import numpy
#from numpy  import *
import numpy as np
import matplotlib.pyplot as plt
#import dislin
#import timeit
#import igraph
#print igraph.__version__
#dir(igraph)
import csv

############################################
# LOAD DATA
############################################

data = np.genfromtxt('../BSF_OSF.out', dtype=None, delimiter=",",skip_header=1,names=None)

years   = np.zeros( (252) )
AABW30S = np.zeros( (252) )
NADW30S = np.zeros( (252) )
NADW30N = np.zeros( (252) )
BSFDP   = np.zeros( (252) )
BSFWG   = np.zeros( (252) )

for i in range(0,252):
  years[i]   = data[i][0] 
  AABW30S[i] = data[i][1]
  NADW30S[i] = data[i][2]
  NADW30N[i] = data[i][3]
  BSFDP[i]   = data[i][4]
  BSFWG[i]   = data[i][5]

del(data)

with open('../POP_SOM_index.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
#    next(readCSV, None)                    # skip the headers
#    row_count = sum(1 for row in readCSV)  # fileObject is your csv.reader
    year      = np.zeros( (251) )
    SOM_index = np.zeros( (251) )
    nr_y      = 0
#    readCSV = csv.DictReader(csvfile, delimiter=',')
    for row in readCSV:
         year[nr_y]      = row[0]
         SOM_index[nr_y] = row[1]
         nr_y            = nr_y + 1


#with open('MXL/XMXL_XMASK_m.csv') as csvfile:
with open('MXL/XMXL_WMASK_m.csv') as csvfile:
#with open('MXL/XMXL_SO30_m.csv') as csvfile:
#with open('MXL/HMXL_XMASK_m.csv') as csvfile:
#with open('MXL/HMXL_WMASK_m.csv') as csvfile:
#with open('MXL/HMXL_SO30_m.csv') as csvfile:
#with open('MXL/TMXL_XMASK_m.csv') as csvfile:
#with open('MXL/TMXL_WMASK_m.csv') as csvfile:
#with open('MXL/TMXL_SO30_m.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV, None)                    # skip the headers
#    row_count = sum(1 for row in readCSV)  # fileObject is your csv.reader
    val  = np.zeros( (12) )
    data = np.zeros( (252) )
    nr_y = 0
    nr_m = 0
#    readCSV = csv.DictReader(csvfile, delimiter=',')
    for row in readCSV:
         val[nr_m]  = row[2]
         nr_m       = nr_m + 1
         if( nr_m==12 ):
          data[nr_y]=np.sum(val,axis=0)/12.0
          nr_m = 0
          nr_y = nr_y + 1
#        print(row)
#        print(row[0])
#        print(row[0],row[1],row[2],)
#        print(row['FirstColumn'])  # Access by column header instead of column number

#data1=data[252-47:252]
data1=data[:]
del(data)

if( 1==2 ):
 plt.figure()
# C,    = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
# AA,   = plt.plot(years,np.abs(AABW30S), 'b', linewidth=3.0, label='C')
# NS,   = plt.plot(years,NADW30S, 'r', linewidth=3.0, label='C')
 SI,   = plt.plot(year, (SOM_index-np.mean(SOM_index[125:251]))/np.std(SOM_index[125:251]), 'g', linewidth=3.0, label='SOM index')
 C,    = plt.plot(years,(data1-np.mean(data1[125:252]))/np.std(data1[125:252]), 'k', linewidth=3.0, label='C')
 NS,   = plt.plot(years,(NADW30S-np.mean(NADW30S[125:252]))/np.std(NADW30S[125:252]), 'r', linewidth=3.0, label='NADW')
 AA,   = plt.plot(years,(np.abs(AABW30S)-np.mean(np.abs(AABW30S[125:252])))/np.std(np.abs(AABW30S[125:252])), 'b', linewidth=3.0, label='AABW')
# NN,   = plt.plot(years,(NADW30N-np.mean(NADW30N))/np.std(NADW30N), 'm', linewidth=3.0, label='C')
 plt.legend(handles=[AA, NS, C, SI])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('b) AABW and NADW volume transports at 30S')
# plt.savefig('PE_spectrum.png')
 plt.savefig('AABW_NADW_moc_v6.0.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

print( np.mean(NADW30S[125:252]) )
print( np.std(NADW30S[125:252]) )
print( np.mean(np.abs(AABW30S[125:252])) )
print( np.std(np.abs(AABW30S[125:252])) )
print( np.nanmean(BSFDP[125:252]) )
print( np.nanstd(BSFDP[125:252]) )
print( np.nanmean(BSFWG[125:252]) )
print( np.nanstd(BSFWG[125:252]) )

if( 1==1 ):
 plt.figure()
# C,    = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
# DP,   = plt.plot(years,BSFDP, 'b', linewidth=3.0, label='DP')
# WG,   = plt.plot(years,BSFWG, 'r', linewidth=3.0, label='WG')
 SI,   = plt.plot(year, (SOM_index-np.mean(SOM_index[125:251]))/np.std(SOM_index[125:251]), 'g', linewidth=3.0, label='SOM index')
 C,    = plt.plot(years,(data1-np.mean(data1[125:252]))/np.std(data1[125:252]), 'k', linewidth=3.0, label='C')
 DP,   = plt.plot(years,(BSFDP-np.nanmean(BSFDP[125:252]))/np.nanstd(BSFDP[125:252]), 'r', linewidth=3.0, label='DP')
 WG,   = plt.plot(years,(BSFWG-np.nanmean(BSFWG[125:252]))/np.nanstd(BSFWG[125:252]), 'b', linewidth=3.0, label='WG')
 plt.legend(handles=[DP, WG, C, SI])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('a) Drake Passage and Weddell Gyre volume transports')
# plt.savefig('PE_spectrum.png')
 plt.savefig('DP_WG_voltran_v6.0.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()


