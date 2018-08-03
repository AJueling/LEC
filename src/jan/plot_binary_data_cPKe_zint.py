# coding: utf-8

# ------------------------------------------------------------------------------
# IMPORTING FUNCTIONALITY
# ------------------------------------------------------------------------------
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib
#matplotlib.use('Agg')                        # Generating matplotlib graphs without a running X server
import matplotlib.pyplot as plt

from read_binary import read_binary_3D, read_binary_2D,read_binary_2D_int
from mapping import generate_lats_lons, shift_field, avg_plots, anomaly_plots

# ------------------------------------------------------------------------------
# PREPARATION
# ------------------------------------------------------------------------------
imt = 3600
jmt = 2400
km  =   42

location    = '/projects/0/samoc/jan/Andree/'
LEC_file_c1 = location+'LEC_bin_5_c1'
LEC_file_c2 = location+'LEC_bin_5_c2'
LEC_file_c3 = location+'LEC_bin_5_c3'
LEC_file_c4 = location+'LEC_bin_5_c4'
grid_file   = '/home/dijkbio/andre/LEC/input/grid.3600x2400.fob.da'

# reads the coordinates and also returns a shift parameter
# this is because the Basemap package likes the data in a [180W,180E] format
# while the binary files start somwhere around 60W or so and then wrap around
lats,lons,shift = generate_lats_lons(grid_file)

# LEC_bin file record numbers
# units are SI
nrec_rPm  =   1  # [J/m^3]
nrec_rPe  =  43
nrec_rKm  =  85
nrec_rKe  = 127

nrec_gPm  = 169  # [W/m^2]
nrec_gPe  = 170
nrec_gKm  = 171
nrec_gKe  = 172
nrec_gPmh = 173
nrec_gPms = 174
nrec_gPeh = 175
nrec_gPes = 176

nrec_cPem = 177  # [W/m^3]
nrec_cKem = 219
nrec_cPKm = 261
nrec_cPKe = 303

nrec_TEMP = 345  # [degC]

nrec_PU_x = 387  # [W/m^3]
nrec_PV_y = 429

# Geometry
geometry1_file = '/home/dijkbio/andre/LEC/input/geometry1'
nrec_TAREA =  5
nrec_UAREA =  6
nrec_DZT   =  7
nrec_DZU   = 49

TAREA = read_binary_2D(geometry1_file,imt,jmt,nrec_TAREA)   # [m^2]
DZT   = read_binary_3D(geometry1_file,imt,jmt,km,nrec_DZT)  # [m^3]

depths = np.genfromtxt('/home/dijkbio/andre/LEC/input/in_depths.42.dat',\
                       delimiter='  ')[:,0]/1e2
for i in range(1,len(depths)):
    depths[i] += depths[i-1]
tdepths = depths
for i in range(len(depths)):
    if i==0: tdepths[i] *=0.5
    else: tdepths[i] = (tdepths[i-1]+tdepths[i])/2

# ------------------------------------------------------------------------------
# READING AND PREPARING DATA
# ------------------------------------------------------------------------------

# Fortran original grid coordinates for tripolar, curvilinear grid:
# longitude (i) 1=250E, 800=30W, 1100=0E, 1500=40E
# latitude  (j) 1=78.5S, 518=55S, 677=45S, 866=30S, 1181=0S
# depth     (k) 1=5m, 9=96.9m, 10=112m, 14=198m, 17=318m, 35=4125, 42=5875m
# beware that in the Northern hemisphere, the grid is not strictly oriented
# North-South / East-West anymore

def vertical_integration(filename,nrec,kmin,kmax,DZ):
    """
    returns a shifted 2D field vertically integrated between kmin and kmax 
    kmin/kmax are depth indices in the Fortran sense
    (i.e. starting at 1 up to and including 42)
    DZ=DZT: integration happens with partial bottom T-cell depths
    """
    print( 'vertical integration from {:4.1f}m to {:4.1f}m'.format(tdepths[kmin-1],tdepths[kmax-1]) )
    field = read_binary_3D(filename,imt,jmt,km,nrec)
    field = np.sum(np.multiply(field[:,:,kmin-1:kmax-1],DZ[:,:,kmin-1:kmax-1]),axis=2)
    field_shifted = shift_field(field,shift)
    return field_shifted

# it does not make sense to do vertical integrals of gPm

gPm_c1_field_shifted = vertical_integration(LEC_file_c1,nrec_cPKe,10,25,DZT)
gPm_c2_field_shifted = vertical_integration(LEC_file_c2,nrec_cPKe,10,25,DZT)
gPm_c3_field_shifted = vertical_integration(LEC_file_c3,nrec_cPKe,10,25,DZT)
gPm_c4_field_shifted = vertical_integration(LEC_file_c4,nrec_cPKe,10,25,DZT)
#gPm_field = read_binary_2D(LEC_file_c1,imt,jmt,nrec_gPm)
#gPm_field_shifted = shift_field(gPm_field,shift)

def horizontal_integration(field,AREA,imin,imax,jmin,jmax):
    """
    integrates an already shifte 2D field horizontally (*m^2) and returns a number
    between imin:imax and jmin:jmax in the new Python convention
    i.e. i: 0=-180E, 1800=0E, 3599=179.9E
    """
    print( 'horizontal integration aprox. [{:4.1f}E,{:4.1f}E]x[{:4.1f}N,{:4.1f}N]'.format(\
           lons[imin,0], lons[imax,0],lats[0,jmin],lats[0,jmax]) )
    integral = np.sum( np.multiply( field[imin:imax,jmin:jmax],\
                                     AREA[imin:imax,jmin:jmax] ) )
    return integral

# example
gPm_int = horizontal_integration(gPm_c1_field_shifted,TAREA,0,imt-1,0,866)
print('cPKe SO30: ', gPm_int/1e9, 'GW')
#exit()
# if( 1==2 ):
#  gPm_field = read_binary_2D(LEC_300_file_c1,imt,jmt,2)  
#  gPm_c1_shifted = shift_field(gPm_field,shift)                 
#  gPm_field = read_binary_2D(LEC_300_file_c2,imt,jmt,2)
#  gPm_c2_shifted = shift_field(gPm_field,shift)              
#  gPm_field = read_binary_2D(LEC_300_file_c3,imt,jmt,2)
#  gPm_c3_shifted = shift_field(gPm_field,shift)              
#  gPm_field = read_binary_2D(LEC_300_file_c4,imt,jmt,2)
#  gPm_c4_shifted = shift_field(gPm_field,shift)              
# print(np.max(gPm_c1_shifted))
# exit()
# gPm_c1_shifted=gPm_c1_shifted

# ------------------------------------------------------------------------------
# PLOTTING
# ------------------------------------------------------------------------------

# if you want to plot rectangluar boxes (that are properly projected onto the map)
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

#def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax,c,name):
def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax,c):
    n=40
    xs = [np.linspace(lonmin,lonmax,n),          np.linspace(lonmax,lonmax,n),          np.linspace(lonmax,lonmin,n),          np.linspace(lonmin,lonmin,n)]
    ys = [np.linspace(latmin,latmin,n),          np.linspace(latmin,latmax,n),          np.linspace(latmax,latmax,n),          np.linspace(latmax,latmin,n)]
    xs = [item for sublist in xs for item in sublist]
    ys = [item for sublist in ys for item in sublist]
#    bmap.plot(xs, ys,latlon = True,color=c,lw=3,label=name)
    bmap.plot(xs, ys,latlon = True,color=c,lw=3)


def rectangle_polygon(m,lonmin,lonmax,latmin,latmax):
    n=40
    xs = [np.linspace(lonmin,lonmax,n),          np.linspace(lonmax,lonmax,n),          np.linspace(lonmax,lonmin,n),          np.linspace(lonmin,lonmin,n)]
    ys = [np.linspace(latmin,latmin,n),          np.linspace(latmin,latmax,n),          np.linspace(latmax,latmax,n),          np.linspace(latmax,latmin,n)]
    xs = [item for sublist in xs for item in sublist]
    ys = [item for sublist in ys for item in sublist]
    poly_coords = []
    for i in range(len(xs)):
        poly_coords.append(m(xs[i],ys[i]))
    return poly_coords


if( 1==1 ):
 f   = plt.figure(figsize=(12,5))
# choose your map projection
 m = Basemap(projection='spaeqd',lon_0=180,boundinglat=0,resolution='c',round=True)  # South Pole projection
# m = Basemap(projection='kav7',lon_0=0,resolution='c')  # global projection
 x, y = m(lons,lats)  # transforms the lats/lons into map coordinates
 m.fillcontinents(color='darkgrey')
 m.drawparallels(np.arange(-90.,91.,30.),fontsize=12) # 
 m.drawmeridians(np.arange(-180.,181.,60.),fontsize=12)
# cs = m.pcolormesh(x,y,gPm_c3_shifted-gPm_c1_shifted,vmin=-4e5,vmax=4e5,cmap='RdBu_r')#,vmin=-0.01,vmax=0.01,shading='flat')
# cs = m.pcolormesh(x,y,gPm_c1_field_shifted*10,cmap='RdBu_r',vmin=-2.,vmax=2.,shading='flat')
 cs = m.pcolormesh(x,y,np.log(np.abs((gPm_c1_field_shifted)+1e-20)),cmap='RdBu_r',vmin=-10.,vmax=1.,shading='flat')
# plot_rectangle(m,-35,80,-78,-50,'orange','WGKP region')
 plot_rectangle(m,-35,80,-78,-50,'orange')
 cb = m.colorbar(cs,location='right',pad="1%",format='%.1f')#,ticks=[-.1,0,.1])
# cb.set_label(label='[dW/m^2]', fontsize=16)
 cb.set_label(label='log(abs([W/m^2]))', fontsize=16)
# plt.legend(fontsize=16)
 f.tight_layout()
 plt.savefig('cPKe_c1_int9_24_log.png')
# plt.savefig('try.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==1 ):
 f   = plt.figure(figsize=(12,5))
# choose your map projection
 m = Basemap(projection='spaeqd',lon_0=180,boundinglat=0,resolution='c',round=True)  # South Pole projection
# m = Basemap(projection='kav7',lon_0=0,resolution='c')  # global projection
 x, y = m(lons,lats)  # transforms the lats/lons into map coordinates
 m.fillcontinents(color='darkgrey')
 m.drawparallels(np.arange(-90.,91.,30.),fontsize=12) # 
 m.drawmeridians(np.arange(-180.,181.,60.),fontsize=12)
# cs = m.pcolormesh(x,y,gPm_c3_shifted-gPm_c1_shifted,vmin=-4e5,vmax=4e5,cmap='RdBu_r')#,vmin=-0.01,vmax=0.01,shading='flat')
# cs = m.pcolormesh(x,y,(gPm_c3_field_shifted - gPm_c1_field_shifted)*10,cmap='RdBu_r',vmin=-1.,vmax=1.,shading='flat')
 cs = m.pcolormesh(x,y,np.log(np.abs((gPm_c3_field_shifted - gPm_c1_field_shifted)+1e-20)),cmap='RdBu_r',vmin=-10.,vmax=1.,shading='flat')
# plot_rectangle(m,-35,80,-78,-50,'orange','WGKP region')
 plot_rectangle(m,-35,80,-78,-50,'orange')
 cb = m.colorbar(cs,location='right',pad="1%",format='%.1f')#,ticks=[-.1,0,.1])
# cb.set_label(label='[dW/m^2]', fontsize=16)
 cb.set_label(label='log(abs([W/m^2]))', fontsize=16)
# plt.legend(fontsize=16)
 f.tight_layout()
 plt.savefig('cPKe_c3_c1_int9_24_log.png')
# plt.savefig('try.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==1 ):
 f   = plt.figure(figsize=(12,5))
# choose your map projection
 m = Basemap(projection='spaeqd',lon_0=180,boundinglat=0,resolution='c',round=True)  # South Pole projection
# m = Basemap(projection='kav7',lon_0=0,resolution='c')  # global projection
 x, y = m(lons,lats)  # transforms the lats/lons into map coordinates
 m.fillcontinents(color='darkgrey')
 m.drawparallels(np.arange(-90.,91.,30.),fontsize=12) # 
 m.drawmeridians(np.arange(-180.,181.,60.),fontsize=12)
# cs = m.pcolormesh(x,y,gPm_c3_shifted-gPm_c1_shifted,vmin=-4e5,vmax=4e5,cmap='RdBu_r')#,vmin=-0.01,vmax=0.01,shading='flat')
# cs = m.pcolormesh(x,y,(gPm_c2_field_shifted - gPm_c1_field_shifted)*10,cmap='RdBu_r',vmin=-1.,vmax=1.,shading='flat')
 cs = m.pcolormesh(x,y,np.log(np.abs((gPm_c2_field_shifted - gPm_c1_field_shifted)+1e-20)),cmap='RdBu_r',vmin=-10.,vmax=1.,shading='flat')
# plot_rectangle(m,-35,80,-78,-50,'orange','WGKP region')
 plot_rectangle(m,-35,80,-78,-50,'orange')
 cb = m.colorbar(cs,location='right',pad="1%",format='%.1f')#,ticks=[-.1,0,.1])
# cb.set_label(label='[dW/m^2]', fontsize=16)
 cb.set_label(label='log(abs([W/m^2]))', fontsize=16)
# plt.legend(fontsize=16)
 f.tight_layout()
 plt.savefig('cPKe_c2_c1_int9_24_log.png')
# plt.savefig('try.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

if( 1==1 ):
 f   = plt.figure(figsize=(12,5))
# choose your map projection
 m = Basemap(projection='spaeqd',lon_0=180,boundinglat=0,resolution='c',round=True)  # South Pole projection
# m = Basemap(projection='kav7',lon_0=0,resolution='c')  # global projection
 x, y = m(lons,lats)  # transforms the lats/lons into map coordinates
 m.fillcontinents(color='darkgrey')
 m.drawparallels(np.arange(-90.,91.,30.),fontsize=12) # 
 m.drawmeridians(np.arange(-180.,181.,60.),fontsize=12)
# cs = m.pcolormesh(x,y,gPm_c3_shifted-gPm_c1_shifted,vmin=-4e5,vmax=4e5,cmap='RdBu_r')#,vmin=-0.01,vmax=0.01,shading='flat')
# cs = m.pcolormesh(x,y,(gPm_c4_field_shifted - gPm_c1_field_shifted)*10,cmap='RdBu_r',vmin=-1.,vmax=1.,shading='flat')
 cs = m.pcolormesh(x,y,np.log(np.abs((gPm_c4_field_shifted - gPm_c1_field_shifted)+1e-20)),cmap='RdBu_r',vmin=-10.,vmax=1.,shading='flat')
# plot_rectangle(m,-35,80,-78,-50,'orange','WGKP region')
 plot_rectangle(m,-35,80,-78,-50,'orange')
 cb = m.colorbar(cs,location='right',pad="1%",format='%.1f')#,ticks=[-.1,0,.1])
# cb.set_label(label='[dW/m^2]', fontsize=16)
 cb.set_label(label='log(abs([W/m^2]))', fontsize=16)
# plt.legend(fontsize=16)
 f.tight_layout()
 plt.savefig('cPKe_c4_c1_int9_24_log.png')
# plt.savefig('try.eps', format='eps', dpi=100)
 plt.show()
 plt.clf()

