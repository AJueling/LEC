# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 13:53:07 2017

@author: andre
"""

import numpy as np
from read_binary import read_binary_2D_double
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap

def generate_lats_lons(grid_file):
    """
    genrates lats and lons fields and shifts them so they are increasing 
    (important for plotting with Basemap)
    """
    imt,jmt = 3600,2400
    lats = read_binary_2D_double(grid_file,imt,jmt,1)
    lons = read_binary_2D_double(grid_file,imt,jmt,2)
    
    shift = np.zeros((jmt),dtype='int')
    for j in range(jmt):
        if j<jmt-1:  b = imt-np.argmin(lons[:,j])
        if j==jmt-1: b = 900
        lats[:,j] = 180/np.pi*np.roll(lats[:,j],b)
        lons[:,j] = 180/np.pi*np.roll(lons[:,j],b)
        shift[j]  = b
    lons[imt-1,jmt-1] = 90.
    
    return lats, lons, shift
    
def shift_field(field,shift):
    """
    shifts a 2D (imt,jmt) field
    """
    imt,jmt = 3600,2400
    shifted = np.zeros((imt,jmt))
    for j in range(jmt):
        shifted[:,j]  = np.roll(field[:,j],shift[j])
    return shifted
    
def vert_int(field,DZT,kmin=0,kmax=41):
    """
    vertically integrates a field using partial or full bottom depths
    """
    vert_int_field = np.sum(np.multiply(field[:,:,kmin:kmax],DZT[:,:,kmin:kmax]),axis=2)
    return vert_int_field
    
    
def avg_plots(field,lons,lats,minv,maxv,name):
    """
    plts the (average) 2D  field on a robin and polar projection
    """
    f   = plt.figure(figsize=(17,6))
    gs  = gridspec.GridSpec(1, 2,width_ratios=[5,3])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    m = Basemap(projection='robin',lon_0=0,resolution='l',ax=ax1)
    m.fillcontinents(color='k')
    x, y = m(lons,lats)
    cs = m.pcolormesh(x,y,field,vmin=minv,vmax=maxv,cmap='RdBu_r',shading='flat')

    m2 = Basemap(projection='spaeqd',boundinglat=-10,lon_0=-60,resolution='l',round=True,ax=ax2)
    cs = m2.pcolormesh(lons,lats,field,vmin=minv,vmax=maxv,cmap='RdBu_r',shading='flat',latlon=True)
    m2.fillcontinents(color='k')
    m2.colorbar(cs,location='right',pad="10%")
    
    f.suptitle(name+' mean [W/m^2]', fontsize=20, fontweight='bold')
    f.subplots_adjust(top=0.9)
    f.tight_layout()
    plt.savefig('Figures/avg_map_'+name)
    
def anomaly_plots(field_avg,field_yrly,lons,lats,minv1,maxv1,minv2,maxv2,name,yr):
    """
    plots mean, 'instantaneous' and difference
    """
    f   = plt.figure(figsize=(17,17))
    gs  = gridspec.GridSpec(3, 2,width_ratios=[5,3])
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[1,0])
    ax4 = plt.subplot(gs[1,1])
    ax5 = plt.subplot(gs[2,0])
    ax6 = plt.subplot(gs[2,1])

    m = Basemap(projection='robin',lon_0=0,resolution='l',ax=ax1)
    m.fillcontinents(color='k')
    x, y = m(lons,lats)
    cs = m.pcolormesh(x,y,field_avg,vmin=minv1,vmax=maxv1,cmap='RdBu_r',shading='flat')

    m2 = Basemap(projection='spaeqd',boundinglat=-10,lon_0=-60,resolution='l',round=True,ax=ax2)
    cs = m2.pcolormesh(lons,lats,field_avg,vmin=minv1,vmax=maxv1,cmap='RdBu_r',shading='flat',latlon=True)
    m2.fillcontinents(color='k')
    m2.colorbar(cs,location='right',pad="7%")

    m = Basemap(projection='robin',lon_0=0,resolution='l',ax=ax3)
    m.fillcontinents(color='k')
    cs = m.pcolormesh(x,y,field_yrly,vmin=minv1,vmax=maxv1,cmap='RdBu_r',shading='flat')

    m2 = Basemap(projection='spaeqd',boundinglat=-10,lon_0=-60,resolution='l',round=True,ax=ax4)
    cs = m2.pcolormesh(lons,lats,field_yrly,vmin=minv1,vmax=maxv1,cmap='RdBu_r',shading='flat',latlon=True)
    m2.fillcontinents(color='k')
    m2.colorbar(cs,location='right',pad="7%")

    m = Basemap(projection='robin',lon_0=0,resolution='l',ax=ax5)
    m.fillcontinents(color='k')
    cs = m.pcolormesh(x,y,field_yrly-field_avg,vmin=minv2,vmax=maxv2,cmap='RdYlBu_r',shading='flat')

    m2 = Basemap(projection='spaeqd',boundinglat=-10,lon_0=-60,resolution='l',round=True,ax=ax6)
    cs = m2.pcolormesh(lons,lats,field_yrly-field_avg,vmin=minv2,vmax=maxv2,cmap='RdYlBu_r',shading='flat',latlon=True)
    m2.fillcontinents(color='k')
    m2.colorbar(cs,location='right',pad="7%")

    f.suptitle(name+' year %3i [W/m^2]'%yr, fontsize=20, fontweight='bold')
    plt.tight_layout()
    f.subplots_adjust(top=0.96)
    plt.savefig('Figures/anomaly_map_'+name+'_%3i'%yr)
    plt.close()
