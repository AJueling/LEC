import numpy as np
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.collections import PatchCollection

from read_output import LEC_global, analyze_LEC

def LEC4_overview(df, name, letter):
    """
    
    """    
    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_aspect('equal')
    ax.set_ylim((0,5))
    ax.set_xlim((0,5))
    grid = np.mgrid[0.5:4.5:5j, 0.5:4.5:5j].reshape(2, -1).T
    patches1, patches2 = [], []

    pos     = [8, 6, 18, 16, 3, 1, 23, 21, 7, 17, 13, 11, 9, 5, 19, 15]
    titles  = ['mean\npotential\nenergy','eddy\npotential\nenergy','mean\nkinetic\nenergy','eddy\nkinetic\nenergy',\
               "G(P$_m$)","G(P$_e$)","G(K$_m$)","G(K$_e$)",\
               "C(P$_e$,P$_m$)","C(K$_e$,K$_m$)","C(P$_m$,K$_m$)","C(P$_e$,K$_e$)",\
               "D(P$_m$)","D(P$_e$)","D(K$_m$)","D(K$_e$)"]
    POP     = [df['rPm'].mean() ,df['rPe'].mean() ,df['rKm'].mean() ,df['rKe'].mean(),\
               df['gPm'].mean() ,df['gPe'].mean() ,df['gKm'].mean() ,df['gKe'].mean(),\
               df['cPem'].mean(),df['cKem'].mean(),df['cPKm'].mean(),df['cPKe'].mean(),\
               df['dPm'].mean() ,df['dPe'].mean() ,df['dKm'].mean() ,df['dKe'].mean()]
    POP_var = [df['rPm'].std()  ,df['rPe'].std()  ,df['rKm'].std()  ,df['rKe'].std(),\
               df['gPm'].std()  ,df['gPe'].std()  ,df['gKm'].std()  ,df['gKe'].std(),\
               df['cPem'].std() ,df['cKem'].std() ,df['cPKm'].std() ,df['cPKe'].std(),\
               df['dPm'].std()  ,df['dPe'].std()  ,df['dKm'].std()  ,df['dKe'].std()]
    total   = df['gPm'].mean() + df['gPe'].mean() + df['gKm'].mean() + df['gKe'].mean()

    color_values2 = []

    rarrow_list = []
    larrow_list = []
    darrow_list = []
    uarrow_list = []

    ax.text(grid[4][0]-.45,grid[4][1]+.25,f'{letter}) {name}',fontsize=24,ha='left')

    
    for i in [4,5,10,11]: #rarrow_list    
        if POP[i]>=0.0:
            rarrow_list.append(i)
        else:
            larrow_list.append(i)
    for i in [6,7]: #rarrow_list    
        if POP[i]>=0.0:
            larrow_list.append(i)
        else:
            rarrow_list.append(i)
    for i in [8,9,12,14]: #uarrow_list    
        if POP[i]>=0.0:
            uarrow_list.append(i)
        else:
            darrow_list.append(i)
    for i in [13,15]: #darrow_list    
        if POP[i]>=0.0:
            darrow_list.append(i)
        else:
            uarrow_list.append(i)
    
    for i in range(len(pos)):
        if i<4:
            rect = mpatches.FancyBboxPatch(grid[pos[i]] - [0.375, 0.375], 0.75, 0.75,\
                               boxstyle=mpatches.BoxStyle("Round", pad=0.1))
        if i in rarrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0]-.47, y=grid[pos[i]][1], dx=.94,dy=0.0,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in larrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0]+.47, y=grid[pos[i]][1], dx=-.94,dy=0.0,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in darrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0], y=grid[pos[i]][1]+.47, dx=0.,dy=-.94,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in uarrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0], y=grid[pos[i]][1]-.47, dx=0.,dy=.94,\
                       width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)

        if i<4:    # energy reservoirs
            plt.text(grid[pos[i]][0], grid[pos[i]][1]+0.2, titles[i],\
                     ha="center",va="center", family='sans-serif', size=16,weight='bold')
            plt.text(grid[pos[i]][0], grid[pos[i]][1]-0.14, "{:4.1f}".format(POP[i]/1e18)+' EJ',\
                     ha="center",va="center", family='sans-serif', size=16)
            plt.text(grid[pos[i]][0], grid[pos[i]][1]-.3, "$\pm${:3.2f}".format(POP_var[i]/1e18)+' EJ',\
                     ha="center",va="center", family='sans-serif', size=16)
            patches1.append(rect)
        elif i>=4: # power transfer terms
            if   i in larrow_list: a,b =  .08, 0.
            elif i in rarrow_list: a,b = -.08, 0.
            elif i in uarrow_list: a,b = 0.  , -.08
            elif i in darrow_list: a,b = 0.  ,  .08
            plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]+0.24+b, titles[i],\
                     ha="center",va="center", family='sans-serif', size=16,weight='bold')
            plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]+0.08+b, "({:4.1f}".format(abs(POP[i]/total)*100.0)+' %)',\
                     ha="center",va="center", family='sans-serif', size=14)
            plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]-0.08+b, "{:4.2f}".format(POP[i]/1e12)+' TW',\
                     ha="center",va="center", family='sans-serif', size=16)
            plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]-0.24+b, "$\pm${:2.0f}".format(abs(POP_var[i])/1e9)+' GW',\
                     ha="center",va="center", family='sans-serif', size=16)
            
            patches2.append(arrw)
            color_values2.append(abs(POP[i]))

    collection1 = PatchCollection(patches1, color='CornflowerBlue', alpha=.6)
    ax.add_collection(collection1)

    collection2 = PatchCollection(patches2, cmap='autumn', alpha=.7)
    collection2.set_array(np.array(-np.array(color_values2)))
    ax.add_collection(collection2)

    plt.axis('off')
    
    plt.savefig('../../results/SOM_paper/LEC4_overview_'+name+'.png',bbox_inches='tight',dpi=100)
    plt.savefig('../../results/SOM_paper/LEC4_overview_'+name+'.eps',bbox_inches='tight',format='eps')



def LEC4_BT_overview(df, name, letter):
    """
    with Km boundary term
    df contains 'df_SO30' and 'budget' dataframes
    """    
    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_aspect('equal')
    ax.set_ylim((0,5))
    ax.set_xlim((0,5))
    grid = np.mgrid[0.5:4.5:5j, 0.5:4.5:5j].reshape(2, -1).T
    patches1, patches2 = [], []

    pos     = [8, 6, 18, 16, 3, 1, 23, 21, 7, 17, 13, 11, 9, 5, 19, 15, 24]
    titles  = ['mean\npotential\nenergy','eddy\npotential\nenergy','mean\nkinetic\nenergy','eddy\nkinetic\nenergy',\
               "G(P$_m$)","G(P$_e$)","G(K$_m$)","G(K$_e$)",\
               "C(P$_e$,P$_m$)","C(K$_e$,K$_m$)","C(P$_m$,K$_m$)","C(P$_e$,K$_e$)",\
               "D/B(P$_m$)","D/B(P$_e$)","D(K$_m$)","D/B(K$_e$)",\
               "B(K$_m$)"]
    POP     = [df['rPm'].mean() ,df['rPe'].mean() ,df['rKm'].mean()     ,df['rKe'].mean(), \
               df['gPm'].mean() ,df['gPe'].mean() ,df['gKm'].mean()     ,df['gKe'].mean(), \
               df['cPem'].mean(),df['cKem'].mean(),df['cPKm'].mean()    ,df['cPKe'].mean(),\
               df['dPm'].mean() ,df['dPe'].mean() ,df['dKm_mbt'].mean() ,df['dKe'].mean(), \
               df['bKm'].mean() ]
    POP_var = [df['rPm'].std()  ,df['rPe'].std()  ,df['rKm'].std()      ,df['rKe'].std(),  \
               df['gPm'].std()  ,df['gPe'].std()  ,df['gKm'].std()      ,df['gKe'].std(),  \
               df['cPem'].std() ,df['cKem'].std() ,df['cPKm'].std()     ,df['cPKe'].std(), \
               df['dPm'].std()  ,df['dPe'].std()  ,df['dKm_mbt'].std()  ,df['dKe'].std(),  \
               df['bKm'].std() ]
    total   = df['gPm'].mean() + df['gPe'].mean() + df['gKm'].mean() + df['gKe'].mean()

    color_values2 = []

    rarrow_list = []
    larrow_list = []
    darrow_list = []
    uarrow_list = []
    barrow_list = [16]

    ax.text(grid[4][0]-.45,grid[4][1]+.25,f'{letter}) {name}',fontsize=24,ha='left')

    
    for i in [4,5,10,11]: #rarrow_list    
        if POP[i]>=0.0:
            rarrow_list.append(i)
        else:
            larrow_list.append(i)
    for i in [6,7]: #rarrow_list    
        if POP[i]>=0.0:
            larrow_list.append(i)
        else:
            rarrow_list.append(i)
    for i in [8,9,12,14]: #uarrow_list    
        if POP[i]>=0.0:
            uarrow_list.append(i)
        else:
            darrow_list.append(i)
    for i in [13,15]: #darrow_list    
        if POP[i]>=0.0:
            darrow_list.append(i)
        else:
            uarrow_list.append(i)
    
    for i in range(len(pos)):
        if i<4:
            rect = mpatches.FancyBboxPatch(grid[pos[i]] - [0.375, 0.375], 0.75, 0.75,\
                               boxstyle=mpatches.BoxStyle("Round", pad=0.1))
        if i in rarrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0]-.47, y=grid[pos[i]][1], dx=.94,dy=0.0,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in larrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0]+.47, y=grid[pos[i]][1], dx=-.94,dy=0.0,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in darrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0], y=grid[pos[i]][1]+.47, dx=0.,dy=-.94,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in uarrow_list:
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0], y=grid[pos[i]][1]-.47, dx=0.,dy=.94,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
        if i in barrow_list: 
            arrw = mpatches.FancyArrow(x=grid[pos[i]][0]-.35, y=grid[pos[i]][1]-.35, dx=0.65,dy=0.65,\
                                   width=0.7,head_width=0.9,head_length=0.2,length_includes_head=True)
            

        if i<4:    # energy reservoirs
            plt.text(grid[pos[i]][0], grid[pos[i]][1]+0.2, titles[i],\
                     ha="center",va="center", family='sans-serif', size=16,weight='bold')
            if POP[i]>1e18: 
                plt.text(grid[pos[i]][0], grid[pos[i]][1]-0.14, "{:4.1f}".format(POP[i]/1e18)+' EJ',\
                        ha="center",va="center", family='sans-serif', size=16)
                plt.text(grid[pos[i]][0], grid[pos[i]][1]-.3, "$\pm${:3.1f}".format(POP_var[i]/1e18)+' EJ',\
                        ha="center",va="center", family='sans-serif', size=16)
                patches1.append(rect)
            elif POP[i]<1e18 and POP[i]>1e15: 
                plt.text(grid[pos[i]][0], grid[pos[i]][1]-0.14, "{:4.1f}".format(POP[i]/1e15)+' PJ',\
                        ha="center",va="center", family='sans-serif', size=16)
                plt.text(grid[pos[i]][0], grid[pos[i]][1]-.3, "$\pm${:3.1f}".format(POP_var[i]/1e15)+' PJ',\
                        ha="center",va="center", family='sans-serif', size=16)
                patches1.append(rect)
            elif POP[i]<1e15 and POP[i]>1e12: 
                plt.text(grid[pos[i]][0], grid[pos[i]][1]-0.14, "{:4.1f}".format(POP[i]/1e12)+' TJ',\
                        ha="center",va="center", family='sans-serif', size=16)
                plt.text(grid[pos[i]][0], grid[pos[i]][1]-.3, "$\pm${:3.1f}".format(POP_var[i]/1e12)+' TJ',\
                        ha="center",va="center", family='sans-serif', size=16)
                patches1.append(rect)
        elif i>=4: 
            if i<16: # power transfer terms
                if   i in larrow_list: a,b =  .08, 0.
                elif i in rarrow_list: a,b = -.08, 0.
                elif i in uarrow_list: a,b = 0.  , -.08
                elif i in darrow_list: a,b = 0.  ,  .08
                plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]+0.24+b, titles[i],\
                         ha="center",va="center", family='sans-serif', size=16,weight='bold')
                plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]+0.08+b, "({:4.1f}".format(abs(POP[i]/total)*100.0)+' %)',\
                         ha="center",va="center", family='sans-serif', size=14)
                plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]-0.08+b, "{:4.0f}".format(POP[i]/1e09)+' GW',\
                         ha="center",va="center", family='sans-serif', size=16)
                if POP_var[i]/1e9>1:
                    plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]-0.24+b, "$\pm${:2.0f}".format(abs(POP_var[i])/1e9)+' GW',\
                            ha="center",va="center", family='sans-serif', size=16)
                else:
                    plt.text(grid[pos[i]][0]+a, grid[pos[i]][1]-0.24+b, "$\pm${:2.1f}".format(abs(POP_var[i])/1e9)+' GW',\
                            ha="center",va="center", family='sans-serif', size=16)
            elif i==16: # boundary term
                a = .16/np.sqrt(2)
                plt.text(grid[pos[i]][0]+a , grid[pos[i]][1]+a, titles[i], rotation=-45,\
                         ha="center",va="center", family='sans-serif', size=16,weight='bold')
                plt.text(grid[pos[i]][0]   , grid[pos[i]][1],\
                         "({:4.1f}".format(abs(POP[i]/total)*100.0)+' %)', rotation=-45,\
                         ha="center",va="center", family='sans-serif', size=14)
                plt.text(grid[pos[i]][0]-a , grid[pos[i]][1]-a,\
                         "{:4.0f}".format(POP[i]/1e9)+' GW', rotation=-45,\
                         ha="center",va="center", family='sans-serif', size=16)
                plt.text(grid[pos[i]][0]-2*a, grid[pos[i]][1]-2*a,\
                         "$\pm${:2.0f}".format(abs(POP_var[i])/1e9)+' GW', rotation=-45,\
                         ha="center",va="center", family='sans-serif', size=16)
                
            patches2.append(arrw)
            color_values2.append(abs(POP[i]))
        
    collection1 = PatchCollection(patches1, color='CornflowerBlue', alpha=.6)
    ax.add_collection(collection1)

    collection2 = PatchCollection(patches2, cmap='autumn', alpha=.7)
    collection2.set_array(np.array(-np.array(color_values2)))
    ax.add_collection(collection2)

    plt.axis('off')
    
    
    plt.savefig('../../results/SOM_paper/LEC4_BT_overview_'+name+'.png',bbox_inches='tight',dpi=100)
    plt.savefig('../../results/SOM_paper/LEC4_BT_overview_'+name+'.eps',bbox_inches='tight',format='eps')
    

start_year=278
end_year = 325
        
for t in range(start_year, end_year):
    fh  = '../../results/analyze_LEC/analysis_LEC_5_'+str(t)+\
    '_SO30.out'
    tmp = pd.read_csv(fh)
    tmp.index = [t]
    if t==start_year: df = tmp.drop([t])
    df = df.append(tmp)
    del tmp, fh


dfg, dfg_anom, dfg_norm = LEC_global('../../results',5)
df_SO30, df_SO30_anom, df_SO30_norm = analyze_LEC('../../results',5,'SO30')
df_WGKP, df_WGKP_anom, df_WGKP_norm = analyze_LEC('../../results',5,'WGKP')


LEC4_overview(dfg,'global','a')
LEC4_BT_overview(df_SO30,'SO30','b')
LEC4_BT_overview(df_WGKP,'WGKP','c')