import csv
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
import cartopy
import cartopy.crs as ccrs
import numpy as np
from netCDF4 import Dataset


""" Fig 1: loading data """

data = np.genfromtxt('../../results/BSF_OSF/BSF_OSF.out', dtype=None, delimiter=",",skip_header=1,names=None)
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

# SOM index
with open('../../results/SOM/POP_SOM_index.csv') as csvfile:
    readCSV   = csv.reader(csvfile, delimiter=',')
    year      = np.zeros( (251) )
    SOM_index = np.zeros( (251) )
    nr_y      = 0
    for row in readCSV:
        year[nr_y]      = row[0]
        SOM_index[nr_y] = row[1]
        nr_y            = nr_y + 1

# convection time series
with open('../../results/MXL/XMXL_WMASK_m.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV, None)
    val  = np.zeros( (12) )
    data = np.zeros( (252) )
    nr_y = 0
    nr_m = 0
    for row in readCSV:
        val[nr_m]  = row[2]
        nr_m       = nr_m + 1
        if( nr_m==12 ):
            data[nr_y]=np.sum(val,axis=0)/12.0
            nr_m = 0
            nr_y = nr_y + 1
data1=data[:]
del(data)


""" Fig 1a: BSF time series """

plt.figure(figsize=(6,4))
# C,    = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
# DP,   = plt.plot(years,BSFDP, 'b', linewidth=3.0, label='DP')
# WG,   = plt.plot(years,BSFWG, 'r', linewidth=3.0, label='WG')
SI,   = plt.plot(year, (SOM_index-np.mean(SOM_index[125:251]))/np.std(SOM_index[125:251]),
                 color=sns.color_palette("tab10")[0], linewidth=2.0, label='SOM index')
C,    = plt.plot(years,(data1-np.mean(data1[125:252]))/np.std(data1[125:252]),
                 'k', linewidth=2.0, label='max MLD')
DP,   = plt.plot(years,(BSFDP-np.nanmean(BSFDP[125:252]))/np.nanstd(BSFDP[125:252]),
                 color=sns.color_palette("tab10")[1], linewidth=2.0, label='Drake Passage')
WG,   = plt.plot(years,(BSFWG-np.nanmean(BSFWG[125:252]))/np.nanstd(BSFWG[125:252]),
                 color=sns.color_palette("tab10")[3], linewidth=2.0, label='Weddell Gyre')
plt.legend(handles=[DP, WG, C, SI])
plt.xlabel('year')
plt.ylabel('normalized amplitude')
plt.title('a) Drake Passage and Weddell Gyre volume transports')
plt.savefig('../../results/SOM_paper/DP_WG_voltran_v8', dpi=100, bbox_inches='tight')
plt.savefig('../../results/SOM_paper/DP_WG_voltran_v8.eps', format='eps', dpi=100, bbox_inches='tight')


""" Fig 1b: OSF timeseries"""

plt.figure(figsize=(6,4))
# C,    = plt.plot(range(278,325),(data1-np.mean(data1))/np.std(data1), 'k', linewidth=3.0, label='C')
# AA,   = plt.plot(years,np.abs(AABW30S), 'b', linewidth=3.0, label='C')
# NS,   = plt.plot(years,NADW30S, 'r', linewidth=3.0, label='C')
SI,   = plt.plot(year, (SOM_index-np.mean(SOM_index[125:251]))/np.std(SOM_index[125:251]),
                 color=sns.color_palette("tab10")[0], linewidth=2.0, label='SOM index')
C,    = plt.plot(years,(data1-np.mean(data1[125:252]))/np.std(data1[125:252]),
                 'k', linewidth=2.0, label='max MLD')
NS,   = plt.plot(years,(NADW30S-np.mean(NADW30S[125:252]))/np.std(NADW30S[125:252]),
                 color=sns.color_palette("tab10")[4], linewidth=2.0, label='NADW')
AA,   = plt.plot(years,(np.abs(AABW30S)-np.mean(np.abs(AABW30S[125:252])))/np.std(np.abs(AABW30S[125:252])),
                 color=sns.color_palette("tab10")[2], linewidth=2.0, label='AABW')
plt.legend(handles=[AA, NS, C, SI])
plt.xlabel('year')
plt.ylabel('normalized amplitude')
plt.title('b) AABW and NADW volume transports at 30$^\circ$S')
plt.savefig('../../results/SOM_paper/AABW_NADW_moc_v8', dpi=100, bbox_inches='tight')
plt.savefig('../../results/SOM_paper/AABW_NADW_moc_v8.eps', format='eps', dpi=100, bbox_inches='tight')




""" Fig 10: loading data """

# use data1 = max MXL time series from Fig. 1
MXL=data1[252-47:252]

#
path = '../../results/analyze_LEC/'
data = np.genfromtxt(path+'analysis_LEC_5_278_SOMS.out', dtype=None, delimiter=",",skip_header=1,names=None)
nr_fields = data.shape[0]
nr_years  = 47
del(data)

data_SO30 = np.zeros( (nr_fields,nr_years) )
data_SO45 = np.zeros( (nr_fields,nr_years) )
data_SO55 = np.zeros( (nr_fields,nr_years) )
data_WGKP = np.zeros( (nr_fields,nr_years) )
data_SOMS = np.zeros( (nr_fields,nr_years) )

for t in range(278,325):
    filename_in   = path+"analysis_LEC_5_%01d_SO30.out" % (t)
    data_SO30[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
    filename_in   = path+"analysis_LEC_5_%01d_SO45.out" % (t)
    data_SO45[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
    filename_in   = path+"analysis_LEC_5_%01d_SO55.out" % (t)
    data_SO55[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
    filename_in   = path+"analysis_LEC_5_%01d_WGKP.out" % (t)
    data_WGKP[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)
    filename_in   = path+"analysis_LEC_5_%01d_SOMS.out" % (t)
    data_SOMS[:,t-278] = np.genfromtxt(filename_in, dtype=None, delimiter=",",skip_header=1,names=None)


def plot_phases_I(ax):
    ax.bar(285.5, 10, 0.5, -1.45, color='lightgrey')
    ax.bar(298.5, 10, 0.5, -1.45, color='lightgrey')
    ax.bar(310.5, 10, 0.5, -1.45, color='lightgrey')
    ax.bar(321.5, 10, 0.5, -1.45, color='lightgrey')
    ax.text(291, 2.1, r'D', fontsize=15)
    ax.text(303, 2.1, r'A', fontsize=15)
    ax.text(315, 2.1, r'B', fontsize=15)
    ax.text(279, 2.1, r'C', fontsize=15)
    ax.text(323, 2.1, r'C', fontsize=15)
def plot_phases_II(ax):
    ax.bar(285.5, 10, 0.5, -2.1, color='lightgrey')
    ax.bar(298.5, 10, 0.5, -2.1, color='lightgrey')
    ax.bar(310.5, 10, 0.5, -2.1, color='lightgrey')
    ax.bar(321.5, 10, 0.5, -2.1, color='lightgrey')
    ax.text(291, 2.1, r'D', fontsize=15)
    ax.text(303, 2.1, r'A', fontsize=15)
    ax.text(315, 2.1, r'B', fontsize=15)
    ax.text(279, 2.1, r'C', fontsize=15)
    ax.text(323, 2.1, r'C', fontsize=15)
ts = range(278,325)


""" Fig 10 a-c: WGKP """

f, ax = plt.subplots(1, 3, figsize=(12,4))

C,   = ax[0].plot(ts,(MXL-np.mean(MXL))/np.std(MXL),
                  'k', linewidth=2.0, label='MLD')
gPm, = ax[0].plot(ts,(np.squeeze(data_WGKP[0,:])-np.mean(np.squeeze(data_WGKP[0,:])))/np.std(np.squeeze(data_WGKP[0,:])),
                  color=sns.color_palette("tab10")[0], linewidth=3.0, label='$G(P_m)$')
gPe, = ax[0].plot(ts,(np.squeeze(data_WGKP[1,:])-np.mean(np.squeeze(data_WGKP[1,:])))/np.std(np.squeeze(data_WGKP[1,:])),
                  color=sns.color_palette("tab10")[1], linewidth=3.0, label='$G(P_e)$')
gKm, = ax[0].plot(ts,(np.squeeze(data_WGKP[2,:])-np.mean(np.squeeze(data_WGKP[2,:])))/np.std(np.squeeze(data_WGKP[2,:])),
                  color=sns.color_palette("tab10")[2], linewidth=3.0, label='$G(K_m)$')
gKe, = ax[0].plot(ts,(np.squeeze(data_WGKP[3,:])-np.mean(np.squeeze(data_WGKP[3,:])))/np.std(np.squeeze(data_WGKP[3,:])),
                  color=sns.color_palette("tab10")[9], linewidth=3.0, label='$G(K_e)$')
plot_phases_I(ax[0])
ax[0].legend(handles=[gPm, gPe, gKm, gKe, C], ncol=3, loc='lower center')
ax[0].set_ylim([-2.4, 2.6])
ax[0].set_xlabel('year', fontsize=12)
ax[0].set_ylabel('normalized amplitude', fontsize=12)
ax[0].set_title('a) Power input in WGKP region', fontsize=12)

C,   = ax[1].plot(ts,(MXL-np.mean(MXL))/np.std(MXL),
                  'k', linewidth=2.0, label='MLD')
rPm, = ax[1].plot(ts,(np.squeeze(data_WGKP[8,:])-np.mean(np.squeeze(data_WGKP[8,:])))/np.std(np.squeeze(data_WGKP[8,:])),
                  color=sns.color_palette("tab10")[0], linewidth=3.0, label='$P_m$')
rPe, = ax[1].plot(ts,(np.squeeze(data_WGKP[9,:])-np.mean(np.squeeze(data_WGKP[9,:])))/np.std(np.squeeze(data_WGKP[9,:])),
                  color=sns.color_palette("tab10")[1], linewidth=3.0, label='$P_e$')
rKm, = ax[1].plot(ts,(np.squeeze(data_WGKP[10,:])-np.mean(np.squeeze(data_WGKP[10,:])))/np.std(np.squeeze(data_WGKP[10,:])),
                  color=sns.color_palette("tab10")[2], linewidth=3.0, label='$K_m$')
rKe, = ax[1].plot(ts,(np.squeeze(data_WGKP[11,:])-np.mean(np.squeeze(data_WGKP[11,:])))/np.std(np.squeeze(data_WGKP[11,:])), 
                  color=sns.color_palette("tab10")[9], linewidth=3.0, label='$K_e$')
plot_phases_I(ax[1])
ax[1].legend(handles=[rPm, rPe, rKm, rKe, C], ncol=3, loc='lower center')
ax[1].set_ylim([-2.4, 2.6])
ax[1].set_xlabel('year', fontsize=12)
ax[1].set_ylabel('normalized amplitude', fontsize=12)
ax[1].set_title('b) Energy reservoirs in WGKP region', fontsize=12)

C,    = ax[2].plot(ts,(MXL-np.mean(MXL))/np.std(MXL),
                   'k', linewidth=2.0, label='MLD')
cPKm, = ax[2].plot(ts,(np.squeeze(data_WGKP[14,:])-np.mean(np.squeeze(data_WGKP[14,:])))/np.std(np.squeeze(data_WGKP[14,:])),
                   color=sns.color_palette("tab10")[0], linewidth=3.0, label='$C(P_m,K_m)$')
cPKe, = ax[2].plot(ts,(np.squeeze(data_WGKP[15,:])-np.mean(np.squeeze(data_WGKP[15,:])))/np.std(np.squeeze(data_WGKP[15,:])),
                   color=sns.color_palette("tab10")[1], linewidth=3.0, label='$C(P_e,K_e)$')
cPem, = ax[2].plot(ts,(np.squeeze(data_WGKP[12,:])-np.mean(np.squeeze(data_WGKP[12,:])))/np.std(np.squeeze(data_WGKP[12,:])),
                   color=sns.color_palette("tab10")[2], linewidth=3.0, label='$C(P_e,P_m)$')
cKem, = ax[2].plot(ts,(np.squeeze(data_WGKP[13,:])-np.mean(np.squeeze(data_WGKP[13,:])))/np.std(np.squeeze(data_WGKP[13,:])),
                   color=sns.color_palette("tab10")[9], linewidth=3.0, label='$C(K_e,K_m)$')

plot_phases_I(ax[2])
ax[2].legend(handles=[cPKm, cPKe, cPem, cKem, C], ncol=3, loc='lower center', handletextpad=0.5, columnspacing=1)
ax[2].set_ylim([-2.4, 2.6])
ax[2].set_xlabel('year', fontsize=12)
ax[2].set_ylabel('normalized amplitude', fontsize=12)
ax[2].set_title('c) Energy conversion in WGKP region', fontsize=12)


plt.tight_layout()

plt.savefig('../../results/SOM_paper/Fig_10_I', dpi=100)
plt.savefig('../../results/SOM_paper/Fig_10_I.eps', format='eps', dpi=100)


""" Fig 10: SO30 """

f, ax = plt.subplots(1, 3, figsize=(12,4))

C,   = ax[0].plot(ts,(MXL-np.mean(MXL))/np.std(MXL),
                  'k', linewidth=1.0, label='MLD')
gPm, = ax[0].plot(ts,(np.squeeze(data_SO30[0,:])-np.mean(np.squeeze(data_SO30[0,:])))/np.std(np.squeeze(data_SO30[0,:])),
                  color=sns.color_palette("tab10")[0],
                  linewidth=3.0,
                  label='$G(P_m)$')
gPe, = ax[0].plot(ts,(np.squeeze(data_SO30[1,:])-np.mean(np.squeeze(data_SO30[1,:])))/np.std(np.squeeze(data_SO30[1,:])),
                  color=sns.color_palette("tab10")[1],
                  linewidth=3.0,
#                   linestyle=':',
                  label='$G(P_e)$')
gKm, = ax[0].plot(ts,(np.squeeze(data_SO30[2,:])-np.mean(np.squeeze(data_SO30[2,:])))/np.std(np.squeeze(data_SO30[2,:])),
                  color=sns.color_palette("tab10")[2],
                  linewidth=3.0,
#                   linestyle='--',
                  label='$G(K_m)$')
gKe, = ax[0].plot(ts,(np.squeeze(data_SO30[3,:])-np.mean(np.squeeze(data_SO30[3,:])))/np.std(np.squeeze(data_SO30[3,:])),
                  color=sns.color_palette("tab10")[9],
                  linewidth=3.0,
#                   linestyle=(0, (3, 1, 1, 1, 1, 1)),
                  label='$G(K_e)$')
plot_phases_II(ax[0])
ax[0].legend(handles=[gPm, gPe, gKm, gKe, C], ncol=3, loc='lower center')
ax[0].set_ylim([-3.2, 2.6])
ax[0].set_xlabel('year', fontsize=12)
ax[0].set_ylabel('normalized amplitude', fontsize=12)
ax[0].set_title('d) Power input in SO30 region', fontsize=12)

C,   = ax[1].plot(ts,(MXL-np.mean(MXL))/np.std(MXL),
                  'k', linewidth=1.0, label='MLD')
rPm, = ax[1].plot(ts,(np.squeeze(data_SO30[8,:])-np.mean(np.squeeze(data_SO30[8,:])))/np.std(np.squeeze(data_SO30[8,:])),
                  color=sns.color_palette("tab10")[0],
                  linewidth=3.0,
                  label='$P_m$')
rPe, = ax[1].plot(ts,(np.squeeze(data_SO30[9,:])-np.mean(np.squeeze(data_SO30[9,:])))/np.std(np.squeeze(data_SO30[9,:])),
                  color=sns.color_palette("tab10")[1],
                  linewidth=3.0,
#                   linestyle=':',
                  label='$P_e$')
rKm, = ax[1].plot(ts,(np.squeeze(data_SO30[10,:])-np.mean(np.squeeze(data_SO30[10,:])))/np.std(np.squeeze(data_SO30[10,:])),
                  color=sns.color_palette("tab10")[2],
                  linewidth=3.0,
#                   linestyle='--',
                  label='$K_m$')
rKe, = ax[1].plot(ts,(np.squeeze(data_SO30[11,:])-np.mean(np.squeeze(data_SO30[11,:])))/np.std(np.squeeze(data_SO30[11,:])), 
                  color=sns.color_palette("tab10")[9],
                  linewidth=3.0,
#                   linestyle=(0, (3, 1, 1, 1, 1, 1)),
                  label='$K_e$')
plot_phases_II(ax[1])
ax[1].legend(handles=[rPm, rPe, rKm, rKe, C], ncol=3, loc='lower center')
ax[1].set_ylim([-3.2, 2.6])
ax[1].set_xlabel('year', fontsize=12)
ax[1].set_ylabel('normalized amplitude', fontsize=12)
ax[1].set_title('e) Energy reservoirs in SO30 region', fontsize=12)

C,    = ax[2].plot(ts,(MXL-np.mean(MXL))/np.std(MXL),
                   'k', linewidth=1.0, label='MLD')
cPKm, = ax[2].plot(ts,(np.squeeze(data_SO30[14,:])-np.mean(np.squeeze(data_SO30[14,:])))/np.std(np.squeeze(data_SO30[14,:])),
                   color=sns.color_palette("tab10")[0],
                   linewidth=3.0,
                   label='$C(P_m,K_m)$')
cPKe, = ax[2].plot(ts,(np.squeeze(data_SO30[15,:])-np.mean(np.squeeze(data_SO30[15,:])))/np.std(np.squeeze(data_SO30[15,:])),
                   color=sns.color_palette("tab10")[1],
                   linewidth=3.0,
#                    linestyle=':',
                   label='$C(P_e,K_e)$')
cPem, = ax[2].plot(ts,(np.squeeze(data_SO30[12,:])-np.mean(np.squeeze(data_SO30[12,:])))/np.std(np.squeeze(data_SO30[12,:])),
                   color=sns.color_palette("tab10")[2],
                   linewidth=3.0,
#                    linestyle='--',
                   label='$C(P_e,P_m)$')
cKem, = ax[2].plot(ts,(np.squeeze(data_SO30[13,:])-np.mean(np.squeeze(data_SO30[13,:])))/np.std(np.squeeze(data_SO30[13,:])),
                   color=sns.color_palette("tab10")[9],
                   linewidth=3.0,
#                    linestyle=(0, (3, 1, 1, 1, 1, 1)),
                   label='$C(K_e,K_m)$')

plot_phases_II(ax[2])
ax[2].legend(handles=[cPKm, cPKe, cPem, cKem, C], ncol=3, loc='lower center', handletextpad=0.5, columnspacing=1)
ax[2].set_ylim([-3.2, 2.6])
ax[2].set_xlabel('year', fontsize=12)
ax[2].set_ylabel('normalized amplitude', fontsize=12)
ax[2].set_title('f) Energy conversion in SO30 region', fontsize=12)


plt.tight_layout()

plt.savefig('../../results/SOM_paper/Fig_10_II', dpi=100)
plt.savefig('../../results/SOM_paper/Fig_10_II.eps', format='eps', dpi=100)


""" Suppl. Figure 2 """

from read_output import LEC_global, analyze_LEC

df_SO30, df_SO30_anom, df_SO30_norm = analyze_LEC('../../results',5,'SO30')
df = df_SO30_anom

def plot_phases_III(ax):
    ax.text(291, 2.4, r'D', fontsize=15)
    ax.text(303, 2.4, r'A', fontsize=15)
    ax.text(315, 2.4, r'B', fontsize=15)
    ax.text(279, 2.4, r'C', fontsize=15)
    ax.text(323, 2.4, r'C', fontsize=15)

f, ax = plt.subplots(1, 2, figsize=(10,4))

gKm,  = ax[0].plot(df.index,  df['gKm'] /1e11            , linewidth=3,
                   color=sns.color_palette("tab10")[0], label='$G(K_m)$')
cKem, = ax[0].plot(df.index,  df['cKem']/1e11            , linewidth=3,
                   color=sns.color_palette("tab10")[1], label='$C(K_m, K_e)$')
cPKm, = ax[0].plot(df.index,  df['cPKm']/1e11            , linewidth=3,
                   color=sns.color_palette("tab10")[2], label='$C(P_m, K_m)$')
dKm,  = ax[0].plot(df.index, -(df['dKm']-df['bKm']) /1e11, linewidth=3,
                   color=sns.color_palette("tab10")[4], label='$-D(K_m)$')  # the dKm in the dataframe still contains the boundary pressure work term
bKm,  = ax[0].plot(df.index, -df['bKm'] /1e11            , linewidth=3,
                   color=sns.color_palette("tab10")[9], label='$-B(K_m)$')

ax[0].set_ylim((-3.5,3))
ax[0].set_xlabel('year', fontsize=12)
ax[0].set_ylabel('power [$10^{11}$ W]')
ax[0].set_title('a) $K_m$ power terms in SO30', fontsize=12)
ax[0].legend(handles=[gKm, cKem, cPKm, dKm, bKm], ncol=3, loc='lower center')
ax[0].bar(285.5, 10, 0.5, -2.2, color='lightgrey')
ax[0].bar(298.5, 10, 0.5, -2.2, color='lightgrey')
ax[0].bar(310.5, 10, 0.5, -2.2, color='lightgrey')
ax[0].bar(321.5, 10, 0.5, -2.2, color='lightgrey')

plot_phases_III(ax[0])

gPm,  = ax[1].plot(df.index,  df['gPm'] /1e11, linewidth=3,
                   color=sns.color_palette("tab10")[0], label='$G(P_m)$')
cPem, = ax[1].plot(df.index,  df['cPem']/1e11, linewidth=3,
                   color=sns.color_palette("tab10")[1], label='$C(P_m, P_e)$')
cPKm, = ax[1].plot(df.index, -df['cPKm']/1e11, linewidth=3,
                   color=sns.color_palette("tab10")[2], label='$-C(P_m, K_m)$')
dPm,  = ax[1].plot(df.index, -df['dPm'] /1e11, linewidth=3,
                   color=sns.color_palette("tab10")[4], label='$-D(P_m)$')

ax[1].set_ylim((-3.5,3))
ax[1].set_xlabel('year', fontsize=12)
ax[1].set_ylabel('power [$10^{11}$ W]')
ax[1].set_title('b) $K_m$ power terms in SO30', fontsize=12)
ax[1].legend(handles=[gPm, cPem, cPKm, dPm], ncol=2, loc='lower right')#, columnspacing=7)
ax[1].bar(285.5, 10, 0.5, -3.5, color='lightgrey')
ax[1].bar(298.5, 10, 0.5, -2.2, color='lightgrey')
ax[1].bar(310.5, 10, 0.5, -2.2, color='lightgrey')
ax[1].bar(321.5, 10, 0.5, -2.2, color='lightgrey')
plot_phases_III(ax[1])


plt.tight_layout()
plt.savefig('../../results/SOM_paper/Suppl_Figure_2')
plt.savefig('../../results/SOM_paper/Suppl_Figure_2.eps', format='eps')



""" Suppl. Figure 3 """

import string
letter_list = list(string.ascii_lowercase)
LEC_term = ['generation', 'reservoir', 'conversion', 'dissipation / boundary']

df_WGKP, df_WGKP_anom, df_WGKP_norm = analyze_LEC('../../results',5,'WGKP')

def plot_phases_IV(ax, loc):
    ax.text(291, loc, r'D', fontsize=14)
    ax.text(303, loc, r'A', fontsize=14)
    ax.text(315, loc, r'B', fontsize=14)
    ax.text(279, loc, r'C', fontsize=14)
    ax.text(323, loc, r'C', fontsize=14)

for i, df_region in enumerate([[df_SO30, df_SO30_anom],[df_WGKP, df_WGKP_anom]]):
    f, ax = plt.subplots(2,4, figsize=(16,8))
    
    for j, df in enumerate(df_region):
        
        
        
        ax[j,0].plot(df.index, df['gPm'], linewidth=3, label='$G(P_m)$')
        ax[j,0].plot(df.index, df['gPe'], linewidth=3, label='$G(P_e)$')
        ax[j,0].plot(df.index, df['gKm'], linewidth=3, label='$G(K_m)$')
        ax[j,0].plot(df.index, df['gKe'], linewidth=3, label='$G(K_e)$')        

        ax[j,1].plot(df.index, df['rPm'], linewidth=3, label='$P_m$')
        ax[j,1].plot(df.index, df['rPe'], linewidth=3, label='$P_e$')
        ax[j,1].plot(df.index, df['rKm'], linewidth=3, label='$K_m$')
        ax[j,1].plot(df.index, df['rKe'], linewidth=3, label='$K_e$')

        ax[j,2].plot(df.index, df['cPKm'], linewidth=3, label='$C(P_m, K_m)$')
        ax[j,2].plot(df.index, df['cPKe'], linewidth=3, label='$C(P_e, K_e)$')
        ax[j,2].plot(df.index, df['cPem'], linewidth=3, label='$C(P_m, P_e)$')
        ax[j,2].plot(df.index, df['cKem'], linewidth=3, label='$C(K_m, K_e)$')

        ax[j,3].plot(df.index, df['dPm']          , linewidth=3, label='$D(P_m)$')
        ax[j,3].plot(df.index, df['dPe']          , linewidth=3, label='$D(P_e)$')
        ax[j,3].plot(df.index, df['dKm']-df['bKm'], linewidth=3, label='$D(K_m)$')
        ax[j,3].plot(df.index, df['dKe']          , linewidth=3, label='$D(K_e)$')
        ax[j,3].plot(df.index, df['bKm']          , linewidth=3, label='$B(K_m)$')
        
        for k in range(4):
            if j==0:
                ax[j,k].plot((278, 324),(0,0), c='k', lw=.5)
            
            letter = letter_list[8*i + 4*j + k]
            ax[j,k].set_title(f'{letter}) {LEC_term[k]}', fontsize=16)
            
            
            ymin, ymax = ax[j,k].get_ylim()
            yextent = ymax - ymin
            ymin = ymax - 1.3*yextent
            ax[j,k].set_ylim((ymin, ymin+1.4*yextent))
            
            plot_phases_IV(ax[j,k], ymax)
            
            ax[j,k].bar(285.5, 1.4/1.3*yextent, 0.5, ymax-yextent, color='lightgrey')
            ax[j,k].bar(298.5, 1.4/1.3*yextent, 0.5, ymax-yextent, color='lightgrey')
            ax[j,k].bar(310.5, 1.4/1.3*yextent, 0.5, ymax-yextent, color='lightgrey')
            ax[j,k].bar(321.5, 1.4/1.3*yextent, 0.5, ymax-yextent, color='lightgrey')
            
            ax[j,k].tick_params(labelsize=15)
            ax[j,k].set_xlabel('year', fontsize=15)
            ax[j,k].legend(ncol=[2,2,2,3][k], fontsize=14, columnspacing=[1,1,.4,.3][k], handlelength=[2,2,1.2,.7][k], loc='lower center')
            if k in [0,2,3]:
                ax[j,k].set_ylabel('power [W]', fontsize=15)
            else:
                ax[j,k].set_ylabel('energy [J]', fontsize=15)

    plt.tight_layout()
    region = ['SO30','WGKP'][i]
    plt.savefig(f'../../results/SOM_paper/Suppl_Figure_3_{region}')
    plt.savefig(f'../../results/SOM_paper/Suppl_Figure_3_{region}.eps', format='eps')