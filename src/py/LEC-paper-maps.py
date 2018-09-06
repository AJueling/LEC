"""
"Energetics of the Southern Ocean Mode" paper submitted to JGR: Oceans
this script creates the maps
it needs the 'nc/LEC_SOM_paper_maps_A_D.nc' that is created with the 'creating_netCDF.ipynb'
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import cartopy
import cartopy.crs as ccrs
import numpy as np
from netCDF4 import Dataset
matplotlib.use('Agg')
#%matplotlib inline


from read_binary import read_binary_2D, read_binary_2D_double, read_binary_3D
from grid import generate_lats_lons, shift_field


def rectangle_polygon(lonmin,lonmax,latmin,latmax):
    n=50
    xs = [np.linspace(lonmin,lonmax,n), np.linspace(lonmax,lonmax,n), np.linspace(lonmax,lonmin,n), np.linspace(lonmin,lonmin,n)]
    ys = [np.linspace(latmin,latmin,n), np.linspace(latmin,latmax,n), np.linspace(latmax,latmax,n), np.linspace(latmax,latmin,n)]
    xs = [item for sublist in xs for item in sublist]
    ys = [item for sublist in ys for item in sublist]
    poly_coords = np.swapaxes(np.array([xs, ys]),0,1)
    return poly_coords


def map_geography(ax):
    """adds land, WGKP ploygon, lat/lon labels"""
    ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='lightgray')
    ax.add_patch(mpatches.Polygon(xy=rectangle_polygon(-35,80,-78,-50),
                                  facecolor='none',
                                  edgecolor='orange',
                                  linewidth=2,
                                  transform=ccrs.PlateCarree()))
    gl = ax.gridlines(xlocs=xticks, ylocs=yticks)
    gl.n_steps = 500
    
    ax.text(180,-62,'60S',
            rotation=5,
            horizontalalignment='left',
            color='darkgray',
            transform=ccrs.PlateCarree())
    ax.text(180,-34,'30S',
            rotation=5,
            horizontalalignment='left',
            color='darkgray',
            transform=ccrs.PlateCarree())

    ax.text(62,-15,'60E',
            rotation=30, horizontalalignment='right', verticalalignment='bottom',
            color='darkgray',
            transform=ccrs.PlateCarree())
    ax.text(119,-10,'120E',
            rotation=-30, horizontalalignment='right', verticalalignment='bottom',
            color='darkgray',
            transform=ccrs.PlateCarree())

    ax.text(-62,-15,'60W',
            rotation=-30, horizontalalignment='left', verticalalignment='bottom',
            color='darkgray',
            transform=ccrs.PlateCarree())
    ax.text(-119,-10,'120W',
            rotation=30, horizontalalignment='left', verticalalignment='bottom',
            color='darkgray',
            transform=ccrs.PlateCarree())

    ax.text(0,-20,'0$^\circ$',
            rotation=90, horizontalalignment='right', verticalalignment='bottom',
            color='darkgray',
            transform=ccrs.PlateCarree())
    ax.text(-180,-10,'180$^\circ$',
            rotation=90, horizontalalignment='right', verticalalignment='bottom',
            color='darkgray',
            transform=ccrs.PlateCarree())

    
""" Geometry """

imt,jmt,km        = 3600, 2400, 42
grid_file       = '../../input/grid.3600x2400.fob.da'
lats,lons,shift   = generate_lats_lons(grid_file)

nrec_DZT = 7
DZT_file = '../../input/geometry1'
DZT = read_binary_3D(DZT_file, imt, jmt, km, nrec_DZT)

thickness = np.genfromtxt('../../input/in_depths.42.dat', delimiter='  ')[:, 0]/1e2
depths    = [thickness[0]]
print(' k     depths thickness')
for i in range(len(thickness)):
    if i>0: depths.append(depths[-1]+thickness[i-1])
    print(f'{i:2d} {depths[i]:10.1f} {thickness[i]:9.1f}')



# lats/lons are redefined below, then being the tranposed of the ones above


""" Figure 2: map of XMXL and WGKP region """


# XMXL_mean_file    = '/projects/0/samoc/jan/Andree/MXL/XMXL_mean'

# XMXL_mean         = read_binary_2D(XMXL_mean_file,imt,jmt,1)
# XMXL_mean_shifted = shift_field(XMXL_mean,shift)

# proj = ccrs.LambertConformal(central_latitude=-60., central_longitude=30, standard_parallels=(-75,-45))

# fig = plt.figure(figsize=(12,10))
# ax = plt.subplot(2,1,1, projection=proj)
# ax.set_extent([-22,72,-35,-85])
# im = ax.pcolormesh(lons.T, lats.T, XMXL_mean_shifted.T, transform=ccrs.PlateCarree())

# ax.add_patch(mpatches.Polygon(xy=rectangle_polygon(-35,80,-78,-50),
#                               facecolor='none',
#                               edgecolor='orange',
#                               linewidth=2,
#                               transform=ccrs.PlateCarree()))

# ax.add_feature(cartopy.feature.LAND, zorder=1, facecolor='lightgray')

# # latitudes
# ax.text(0   , .43,  '45$^\circ$S', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
# ax.text(.1  , 0  ,  '60$^\circ$S', horizontalalignment='left' , verticalalignment='top'   , transform=ax.transAxes)
# ax.text(.3  , 0  ,  '75$^\circ$S', horizontalalignment='left' , verticalalignment='top'   , transform=ax.transAxes)
# ax.text(.74 , 0  ,  '75$^\circ$S', horizontalalignment='left' , verticalalignment='top'   , transform=ax.transAxes)
# ax.text(.94 , 0  ,  '60$^\circ$S', horizontalalignment='left' , verticalalignment='top'   , transform=ax.transAxes)
# ax.text(1   , .64,  '45$^\circ$S', horizontalalignment='left' , verticalalignment='bottom', transform=ax.transAxes)
# 
# # longitutes
# ax.text(0   , .07,  '60$^\circ$W', horizontalalignment='right', transform=ax.transAxes)
# ax.text(0   , .63,  '30$^\circ$W', horizontalalignment='right', transform=ax.transAxes)
# ax.text(.23 , 1  ,   '0$^\circ$E', horizontalalignment='left' , transform=ax.transAxes)
# ax.text(.53 , 1  ,  '30$^\circ$E', horizontalalignment='left' , transform=ax.transAxes)
# ax.text(.835, 1  ,  '60$^\circ$E', horizontalalignment='left' , transform=ax.transAxes)
# ax.text(1   , .51,  '90$^\circ$E', horizontalalignment='left' , transform=ax.transAxes)
# ax.text(1   , .04, '120$^\circ$E', horizontalalignment='left' , transform=ax.transAxes)

# xticks = np.linspace(-180,180,13,endpoint=True)
# yticks = np.linspace(-90,90,13)
# gl = ax.gridlines(xlocs=xticks, ylocs=yticks, color='k')
# gl.n_steps = 500

# ax2 = plt.subplot(2,1,2)

# # all temperature in WGKP is incomplete, could be recalculated
# # for now I just edited the existing plots into one manually


proj   = ccrs.Orthographic(central_latitude=-90)
xticks = np.linspace(-180,180,7,endpoint=True)
yticks = np.linspace(-90,90,7)

""" Bathymetry map """    
bathymetry = np.sum(DZT, axis=2)
bathymetry = shift_field(bathymetry, shift)

fig = plt.figure(figsize=(7, 8))
ax = plt.axes(projection=proj)
cax, kw = matplotlib.colorbar.make_axes(ax, location='bottom', pad=0.05, shrink=.8)
im = ax.pcolormesh(lons.T, lats.T, bathymetry.T/1000, transform=ccrs.PlateCarree(), cmap='YlGnBu', vmax=6)#, norm=norm)
map_geography(ax)
cbar = fig.colorbar(im, cax=cax,**kw)
cbar.ax.tick_params(labelsize=14)
cbar.set_ticks(np.linspace(0,6,7))
label = cbar.set_label('bathymetry [km]', size=16)
plt.savefig('../../results/SOM_paper/bathymetry', bbox_inches='tight')



""" Figure 4-9: LEC spatial plots """



netCDF_file = '../../results/SOM_paper/nc/LEC_SOM_paper_maps_A_D.nc'
rootgrp = Dataset(netCDF_file, "r", format="NETCDF4")

lons = rootgrp['lon'][:,:]
lats = rootgrp['lat'][:,:]

phases = ['phase_D/', 'phase_A/', 'phase_B/', 'phase_C/',]

phase_label = ['reference\nphase D',
               'phase A-D\nanomaly',
               'phase B-D\nanomaly',
               'phase C-D\nanomaly',
              ]

quantities = ['rPm_vint_100_1500','rPe_vint_100_1500','rKm_vint_100_1500','rKe_vint_100_1500',
              'gPm','gPe','gKm','gKe',
              'cPKm_vint_100_1500','cPKe_vint_100_1500','cPem_vint_100_1500','cKem_vint_100_1500',
             ]

# reference phase plots

ref_cmaps = ['viridis','viridis','viridis','viridis',
             'Spectral_r','Spectral_r','Spectral_r','Spectral_r',
             'viridis','viridis','viridis','viridis',
            ]

ref_range = [[10,15],[6,13],[5,12],[6,13],
             [-.5,.5],[-.5,.5],[-.5,.5],[-.5,.5],
             [-10,1],[-10,1],[-10,1],[-10,1],
            ]

ref_labels = ['log([J/m$^2$])','log([J/m$^2$])','log([J/m$^2$])','log([J/m$^2$])',
              '[dW/m$^2$]','[cW/m$^2$]','[dW/m$^2$]','[cW/m$^2$]',
              'log(abs([W/m$^2$]))','log(abs([W/m$^2$]))','log(abs([W/m$^2$]))','log(abs([W/m$^2$]))',
             ]

ref_factors = [1,1,1,1,
               1e1, 1e2, 1e1, 1e2,
               1,1,1,1,
               ]

# phase anomaly plots

anom_range = [.8,10,5,10,
              1,1,1,1,
              [-10,1],[-10,1],[-10,1],[-10,1],
             ]

anom_labels = ['[MJ/m$^2$]','[kJ/m$^2$]','[kJ/m$^2$]','[kJ/m$^2$]',
               '[cW/m$^2$]','[mW/m$^2$]','[cW/m$^2$]','[mW/m$^2$]',
               'log(abs([W/m$^2$]))','log(abs([W/m$^2$]))','log(abs([W/m$^2$]))','log(abs([W/m$^2$]))',
              ]

anom_factors = [1e-6, 1e-3, 1e-3, 1e-3,
                100,1000,100,1000,
                1,1,1,1,
               ]

letters = [['a)', 'b)', 'c)', 'd)'],
           ['e)', 'f)', 'g)', 'h)']
          ]



for j, quant in enumerate(quantities):  # all quantities
    for i, phase in enumerate(phases):  # phases D, A, B, C
        print(f'j={j} i={i}')
        
        fig = plt.figure(figsize=(7, 8))
        ax  = plt.axes(projection=proj)
        cax, kw = matplotlib.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.8)
        
        if i==0:  # reference phase
            ref_field = rootgrp[phase+quant][:,:]
            the_norm = matplotlib.colors.Normalize(vmin=ref_range[j][0],vmax=ref_range[j][1])
            if j in range(4):  # reservoirs: log
                im = ax.pcolormesh(lons, lats, np.log(ref_field),
                                   norm=the_norm,
                                   cmap=ref_cmaps[j],
                                   transform=ccrs.PlateCarree()
                                  )
            elif j in range(4,8):  # generation: actual fields
                im = ax.pcolormesh(lons, lats, ref_field*ref_factors[j],
                                   norm=the_norm,
                                   cmap=ref_cmaps[j],
                                   transform=ccrs.PlateCarree()
                                  )
            elif j in range(8,12):  # exchange: log(abs(...))
                im = ax.pcolormesh(lons, lats, np.log(abs(ref_field)),
                                   norm=the_norm,
                                   cmap=ref_cmaps[j],
                                   transform=ccrs.PlateCarree()
                                  )
            cbar = fig.colorbar(im,cax=cax,extend='both',**kw)
            cbar.ax.tick_params(labelsize=14)
            label = cbar.set_label(ref_labels[j], size=16)

        else:  # anomaly phases
            anom_field = rootgrp[phase+quant][:,:] - ref_field
            if j<8:  # reservoirs/generation true anomalies
                the_norm = matplotlib.colors.Normalize(vmin=-anom_range[j],vmax=anom_range[j])
                im = ax.pcolormesh(lons, lats, anom_field*anom_factors[j],
                                   norm=the_norm,
                                   cmap='RdBu_r',
                                   transform=ccrs.PlateCarree()
                                  )
                cbar = fig.colorbar(im,cax=cax,extend='both',**kw)
                cbar.ax.tick_params(labelsize=14)
                cbar.set_ticks(np.linspace(-anom_range[j], anom_range[j], 5))
            else:  # conversion: log(abs())
                the_norm = matplotlib.colors.Normalize(vmin=anom_range[j][0],vmax=anom_range[j][1])
                im = ax.pcolormesh(lons, lats, np.log(abs(anom_field)),
                                   norm=the_norm,
                                   cmap='magma',
                                   transform=ccrs.PlateCarree()
                                  )
                cbar = fig.colorbar(im,cax=cax,extend='both',**kw)
                cbar.ax.tick_params(labelsize=14)
            label = cbar.set_label(anom_labels[j], size=16)

        map_geography(ax)

        ax.text(.03, .93, letters[j%2][i], fontsize=20, horizontalalignment='left', transform=ax.transAxes)
        ax.text(.98, .93, phase_label[i], fontsize=16, horizontalalignment='right', transform=ax.transAxes)
        
        #plt.savefig('../../results/SOM_paper/tif/'+quant+'_'+['D','A-D','B-D','C-D'][i]+'.tif', format='tiff', dpi=300, bbox_inches='tight')
        plt.savefig('../../results/SOM_paper/'+quant+'_'+['D','A-D','B-D','C-D'][i], dpi=250, bbox_inches='tight')
        #plt.savefig('../../results/SOM_paper/'+quant+'_'+['D','A-D','B-D','C-D'][i]+'.eps', format='eps') #, bbox_inches='tight')
        #plt.savefig('../../results/SOM_paper/'+quant+'_'+['D','A-D','B-D','C-D'][i]+'.pdf', format='pdf') #, bbox_inches='tight')
        plt.close()