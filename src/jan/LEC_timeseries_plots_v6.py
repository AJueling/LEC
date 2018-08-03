
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

data1=data[252-47:252]
del(data)


data = np.genfromtxt('LEC_results_analyze_LEC/analysis_LEC_5_278_SOMS.out', dtype=None, delimiter=",",skip_header=1,names=None)

nr_fields = data.shape[0]
nr_years  = 47
del(data)

data_SO30 = np.zeros( (nr_fields,nr_years) )
data_SO45 = np.zeros( (nr_fields,nr_years) )
data_SO55 = np.zeros( (nr_fields,nr_years) )
data_WGKP = np.zeros( (nr_fields,nr_years) )
data_SOMS = np.zeros( (nr_fields,nr_years) )

for t in range(278,325):
 filename_in   = "LEC_results_analyze_LEC/analysis_LEC_5_%01d_SO30.out" % (t)
 data_SO30[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
 filename_in   = "LEC_results_analyze_LEC/analysis_LEC_5_%01d_SO45.out" % (t)
 data_SO45[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
 filename_in   = "LEC_results_analyze_LEC/analysis_LEC_5_%01d_SO55.out" % (t)
 data_SO55[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
 filename_in   = "LEC_results_analyze_LEC/analysis_LEC_5_%01d_WGKP.out" % (t)
 data_WGKP[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
 filename_in   = "LEC_results_analyze_LEC/analysis_LEC_5_%01d_SOMS.out" % (t)
 data_SOMS[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)

if( 1==2 ):
 plt.figure()
 C,   = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
 gKe, = plt.plot(range(278,325),(np.squeeze(data_WGKP[3,:])-np.mean(np.squeeze(data_WGKP[3,:])))/np.std(np.squeeze(data_WGKP[3,:])), 'm', linewidth=3.0, label='gKe')
 gPe, = plt.plot(range(278,325),(np.squeeze(data_WGKP[1,:])-np.mean(np.squeeze(data_WGKP[1,:])))/np.std(np.squeeze(data_WGKP[1,:])), 'c', linewidth=3.0, label='gPe')
 gKm, = plt.plot(range(278,325),(np.squeeze(data_WGKP[2,:])-np.mean(np.squeeze(data_WGKP[2,:])))/np.std(np.squeeze(data_WGKP[2,:])), 'r', linewidth=3.0, label='gKm')
 gPm, = plt.plot(range(278,325),(np.squeeze(data_WGKP[0,:])-np.mean(np.squeeze(data_WGKP[0,:])))/np.std(np.squeeze(data_WGKP[0,:])), 'b', linewidth=3.0, label='gPm')
 plt.bar(285.5, 5, 0.5, -2, color='grey')
 plt.bar(300.5, 5, 0.5, -2, color='grey')
 plt.bar(310.5, 5, 0.5, -2, color='grey')
 plt.bar(320.5, 5, 0.5, -2, color='grey')
 plt.text(292, 2.3, r'p1', fontsize=15)
 plt.text(304, 2.3, r'p2', fontsize=15)
 plt.text(314, 2.3, r'p3', fontsize=15)
 plt.text(279, 2.3, r'p4', fontsize=15)
 plt.legend(handles=[C, gPm, gPe, gKm, gKe],loc='lower right')
 plt.ylim([-1.8, 2.6])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('a) Power input in WGKP region')
# plt.savefig('PE_spectrum.png')
 plt.savefig('/Users/janviebahn/Desktop/new/2018/03_job/02_imau/01_SOM_paper/Henk/Jan/pics/input_WGKP_v6.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==2 ):
 plt.figure()
 C,   = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
 rKe, = plt.plot(range(278,325),(np.squeeze(data_WGKP[11,:])-np.mean(np.squeeze(data_WGKP[11,:])))/np.std(np.squeeze(data_WGKP[11,:])), 'm', linewidth=3.0, label='rKe')
 rPe, = plt.plot(range(278,325),(np.squeeze(data_WGKP[9,:])-np.mean(np.squeeze(data_WGKP[9,:])))/np.std(np.squeeze(data_WGKP[9,:])), 'c', linewidth=3.0, label='rPe')
 rKm, = plt.plot(range(278,325),(np.squeeze(data_WGKP[10,:])-np.mean(np.squeeze(data_WGKP[10,:])))/np.std(np.squeeze(data_WGKP[10,:])), 'r', linewidth=3.0, label='rKm')
 rPm, = plt.plot(range(278,325),(np.squeeze(data_WGKP[8,:])-np.mean(np.squeeze(data_WGKP[8,:])))/np.std(np.squeeze(data_WGKP[8,:])), 'b', linewidth=3.0, label='rPm')
 plt.bar(285.5, 5, 0.5, -2, color='grey')
 plt.bar(300.5, 5, 0.5, -2, color='grey')
 plt.bar(310.5, 5, 0.5, -2, color='grey')
 plt.bar(320.5, 5, 0.5, -2, color='grey')
 plt.text(292, 2.3, r'p1', fontsize=15)
 plt.text(304, 2.3, r'p2', fontsize=15)
 plt.text(314, 2.3, r'p3', fontsize=15)
 plt.text(279, 2.3, r'p4', fontsize=15)
 plt.legend(handles=[C, rPm, rPe, rKm, rKe],loc='lower right')
 plt.ylim([-1.8, 2.6])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('b) Energy reservoirs in WGKP region')
# plt.savefig('PE_spectrum.png')
 plt.savefig('/Users/janviebahn/Desktop/new/2018/03_job/02_imau/01_SOM_paper/Henk/Jan/pics/reservoir_WGKP_v6.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==2 ):
 plt.figure()
 C,    = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
 cKem, = plt.plot(range(278,325),(np.squeeze(data_WGKP[13,:])-np.mean(np.squeeze(data_WGKP[13,:])))/np.std(np.squeeze(data_WGKP[13,:])), 'm', linewidth=3.0, label='cKem')
 cPem, = plt.plot(range(278,325),(np.squeeze(data_WGKP[12,:])-np.mean(np.squeeze(data_WGKP[12,:])))/np.std(np.squeeze(data_WGKP[12,:])), 'c', linewidth=3.0, label='cPem')
 cPKe, = plt.plot(range(278,325),(np.squeeze(data_WGKP[15,:])-np.mean(np.squeeze(data_WGKP[15,:])))/np.std(np.squeeze(data_WGKP[15,:])), 'r', linewidth=3.0, label='cPKe')
 cPKm, = plt.plot(range(278,325),(np.squeeze(data_WGKP[14,:])-np.mean(np.squeeze(data_WGKP[14,:])))/np.std(np.squeeze(data_WGKP[14,:])), 'b', linewidth=3.0, label='cPKm')
 plt.bar(285.5, 5, 0.5, -2, color='grey')
 plt.bar(300.5, 5, 0.5, -2, color='grey')
 plt.bar(310.5, 5, 0.5, -2, color='grey')
 plt.bar(320.5, 5, 0.5, -2, color='grey')
 plt.text(292, 2.3, r'p1', fontsize=15)
 plt.text(304, 2.3, r'p2', fontsize=15)
 plt.text(314, 2.3, r'p3', fontsize=15)
 plt.text(279, 2.3, r'p4', fontsize=15)
 plt.legend(handles=[C, cPKm, cPem, cPKe, cKem],loc='lower right')
 plt.ylim([-1.8, 2.6])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('c) Energy conversion in WGKP region')
# plt.savefig('PE_spectrum.png')
 plt.savefig('/Users/janviebahn/Desktop/new/2018/03_job/02_imau/01_SOM_paper/Henk/Jan/pics/conversion_WGKP_v6.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==1 ):
 plt.figure()
 C,   = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
 gKe, = plt.plot(range(278,325),(np.squeeze(data_SO30[3,:])-np.mean(np.squeeze(data_SO30[3,:])))/np.std(np.squeeze(data_SO30[3,:])), 'm', linewidth=3.0, label='gKe')
 gPe, = plt.plot(range(278,325),(np.squeeze(data_SO30[1,:])-np.mean(np.squeeze(data_SO30[1,:])))/np.std(np.squeeze(data_SO30[1,:])), 'c', linewidth=3.0, label='gPe')
 gKm, = plt.plot(range(278,325),(np.squeeze(data_SO30[2,:])-np.mean(np.squeeze(data_SO30[2,:])))/np.std(np.squeeze(data_SO30[2,:])), 'r', linewidth=3.0, label='gKm')
 gPm, = plt.plot(range(278,325),(np.squeeze(data_SO30[0,:])-np.mean(np.squeeze(data_SO30[0,:])))/np.std(np.squeeze(data_SO30[0,:])), 'b', linewidth=3.0, label='gPm')
 plt.bar(285.5, 7, 0.5, -3.5, color='grey')
 plt.bar(300.5, 7, 0.5, -3.5, color='grey')
 plt.bar(310.5, 7, 0.5, -3.5, color='grey')
 plt.bar(320.5, 7, 0.5, -3.5, color='grey')
 plt.text(292, 2.3, r'p1', fontsize=15)
 plt.text(304, 2.3, r'p2', fontsize=15)
 plt.text(314, 2.3, r'p3', fontsize=15)
 plt.text(279, 2.3, r'p4', fontsize=15)
 plt.legend(handles=[C, gPm, gPe, gKm, gKe],loc='lower right')
 plt.ylim([-3.2, 2.7])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('d) Power input in SO30 region')
# plt.savefig('PE_spectrum.png')
 plt.savefig('/Users/janviebahn/Desktop/new/2018/03_job/02_imau/01_SOM_paper/Henk/Jan/pics/input_SO30_v6.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==1 ):
 plt.figure()
 C,   = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
 rKe, = plt.plot(range(278,325),(np.squeeze(data_SO30[11,:])-np.mean(np.squeeze(data_SO30[11,:])))/np.std(np.squeeze(data_SO30[11,:])), 'm', linewidth=3.0, label='rKe')
 rPe, = plt.plot(range(278,325),(np.squeeze(data_SO30[9,:])-np.mean(np.squeeze(data_SO30[9,:])))/np.std(np.squeeze(data_SO30[9,:])), 'c', linewidth=3.0, label='rPe')
 rKm, = plt.plot(range(278,325),(np.squeeze(data_SO30[10,:])-np.mean(np.squeeze(data_SO30[10,:])))/np.std(np.squeeze(data_SO30[10,:])), 'r', linewidth=3.0, label='rKm')
 rPm, = plt.plot(range(278,325),(np.squeeze(data_SO30[8,:])-np.mean(np.squeeze(data_SO30[8,:])))/np.std(np.squeeze(data_SO30[8,:])), 'b', linewidth=3.0, label='rPm')
 plt.bar(285.5, 6, 0.5, -2.5, color='grey')
 plt.bar(300.5, 6, 0.5, -2.5, color='grey')
 plt.bar(310.5, 6, 0.5, -2.5, color='grey')
 plt.bar(320.5, 6, 0.5, -2.5, color='grey')
 plt.text(292, 2.3, r'p1', fontsize=15)
 plt.text(304, 2.3, r'p2', fontsize=15)
 plt.text(314, 2.3, r'p3', fontsize=15)
 plt.text(279, 2.3, r'p4', fontsize=15)
 plt.legend(handles=[C, rPm, rPe, rKm, rKe],loc='lower left')
 plt.ylim([-2.3, 2.6])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('e) Energy reservoirs in SO30 region')
# plt.savefig('PE_spectrum.png')
 plt.savefig('/Users/janviebahn/Desktop/new/2018/03_job/02_imau/01_SOM_paper/Henk/Jan/pics/reservoir_SO30_v6.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==1 ):
 plt.figure()
 C,    = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
 cKem, = plt.plot(range(278,325),(np.squeeze(data_SO30[13,:])-np.mean(np.squeeze(data_SO30[13,:])))/np.std(np.squeeze(data_SO30[13,:])), 'm', linewidth=3.0, label='cKem')
 cPem, = plt.plot(range(278,325),(np.squeeze(data_SO30[12,:])-np.mean(np.squeeze(data_SO30[12,:])))/np.std(np.squeeze(data_SO30[12,:])), 'c', linewidth=3.0, label='cPem')
 cPKe, = plt.plot(range(278,325),(np.squeeze(data_SO30[15,:])-np.mean(np.squeeze(data_SO30[15,:])))/np.std(np.squeeze(data_SO30[15,:])), 'r', linewidth=3.0, label='cPKe')
 cPKm, = plt.plot(range(278,325),(np.squeeze(data_SO30[14,:])-np.mean(np.squeeze(data_SO30[14,:])))/np.std(np.squeeze(data_SO30[14,:])), 'b', linewidth=3.0, label='cPKm')
 plt.bar(285.5, 6, 0.5, -3, color='grey')
 plt.bar(300.5, 6, 0.5, -3, color='grey')
 plt.bar(310.5, 6, 0.5, -3, color='grey')
 plt.bar(320.5, 6, 0.5, -3, color='grey')
 plt.text(292, 2.3, r'p1', fontsize=15)
 plt.text(304, 2.3, r'p2', fontsize=15)
 plt.text(314, 2.3, r'p3', fontsize=15)
 plt.text(279, 2.3, r'p4', fontsize=15)
 plt.legend(handles=[C, cPKm, cPem, cPKe, cKem],loc='lower right')
 plt.ylim([-3, 2.7])
# plt.xscale('log', basex=10)
 plt.xlabel('year')
 plt.ylabel('normalized amplitude')
 plt.title('f) Energy conversion in SO30 region')
# plt.savefig('PE_spectrum.png')
 plt.savefig('/Users/janviebahn/Desktop/new/2018/03_job/02_imau/01_SOM_paper/Henk/Jan/pics/conversion_SO30_v6.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()




