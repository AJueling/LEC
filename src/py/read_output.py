# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 16:32:17 2017

@author: andre
"""
import pandas as pd


def LEC_global(location,ntavg):
    """
    reads gloabl integral output from LEC.f90
    returns pandas dataframe time series
    """
    if ntavg==5:
        start_year=278
        end_year = 325
    
    for t in range(start_year, end_year):
        fh  = location+'/global/global_LEC_'+str(ntavg)+'_%3d.out' % t
        tmp = pd.read_csv(fh)
        tmp.index = [t]
        if t==start_year: dfg = tmp.drop([t])
        dfg = dfg.append(tmp)
        del tmp, fh        

    calculate_total(dfg)
    time_derivative(dfg)
    dfg_anom, dfg_norm = anomalies_norm(dfg)
    return dfg, dfg_anom, dfg_norm

def LEC_levels(location,ntavg):
    """
    reads surface integral output from LEC.f90
    returns pandas dataframe time series
    """
    if ntavg==5:
        start_year=278
        end_year = 325
    
    for t in range(start_year, end_year):
        fh  = location+'/levels/level_LEC_'+str(ntavg)+'_%3d.out' % t
        tmp = pd.read_csv(fh)
        tmp.index = [t]
        if t==start_year: dfl = tmp.drop([t])
        dfl = dfl.append(tmp) 
        del tmp, fh
    
    return dfl


def analyze_LEC(location,ntavg,area_name):
    """
    additionally it calculates the dissipation terms
    returns pandas dataframe time series
    """
    
    if ntavg==5:
        start_year=278
        end_year = 325
        
    for t in range(start_year, end_year):
        fh  = location+'/analyze_LEC/analysis_LEC_'+str(ntavg)+'_'+str(t)+\
        '_'+area_name+'.out'
        tmp = pd.read_csv(fh)
        tmp.index = [t]
        if t==start_year: df = tmp.drop([t])
        df = df.append(tmp)
        del tmp, fh
    
    df['dPm'] = df['gPm'] - df['cPKm'] + df['cPem']
    df['dPe'] = df['gPe'] - df['cPKe'] - df['cPem']
    df['dKm'] = df['gKm'] + df['cPKm'] + df['cKem']
    df['dKe'] = df['gKe'] + df['cPKe'] - df['cKem']
    
    df['bKm'] = df['PU_x_int'] + df['PV_y_int']

    df['dKm_mbt'] = df['dKm'] - df['bKm']  #

    calculate_total(df)
    time_derivative(df)
    df_anom, df_norm = anomalies_norm(df)
    return df, df_anom, df_norm

def calculate_total(df):
    """
    adds the total=mean+eddy to the dataframe
    """    
    quant = ['rP','rK','gP','gK','cPK','dP','dK']
    for q in quant:
        df[q+'t'] = df[q+'m'] + df[q+'e']
    return df
    
def time_derivative(df):
    """
    calculates the time derivative of the energy terms in a dataframe
    assumes equidistant time steps of 1 year
    results in [W]
    """    
    nt = len(df)        
    for E in ['P','K']:
        for s in ['m','e','t']:
            df['r'+E+s+'_dt'] = pd.Series(index=df.index)
            for i in range(nt):
                t = df.index[i]
                if i==0:
                    df1 = df['r'+E+s][t]
                    df2 = df['r'+E+s][t+1]
                elif i==nt-1:
                    df1 = df['r'+E+s][t-1]
                    df2 = df['r'+E+s][t]
                else:
                    df1 = df['r'+E+s][t-1]/2.0
                    df2 = df['r'+E+s][t+1]/2.0
                df['r'+E+s+'_dt'][t] = (df2-df1)/365/24/3600
                del df1, df2

def anomalies_norm(df):
    """
    creates anomaly and normalized [-1,+1] dataframes
    """
    df_anom = df - df.mean()
    df_norm = 2*df / (df.max() - df.min()) - (2*df / (df.max() - df.min())).max() + 1
    return df_anom, df_norm

def SOM_netcdf():
    """
    """
    fh = Dataset('/Users/andre/Downloads/SOM_POP.nc', mode='r')
    # SOM dataframe
    dfs = pd.DataFrame(index=range(1,326))
    
    time = fh.variables['time'][:].squeeze()
    for t in range(len(time)):
        time[t] = int(num2date(time[t],'hours since 0001-01-01 00:00:00').year)
        
    depth     = fh.variables['depth'][:] # in [m]
    tlons     = fh.variables['tlon'][:]
    tlats     = fh.variables['tlat'][:]
    ulons     = fh.variables['ulon'][:]
    ulats     = fh.variables['ulat'][:]
    ni,nj     = len(tlats),len(tlons)
    
    dfs['BSF_drake']     = pd.Series(fh.variables['BSF_drake'],index=time)
    dfs['BSF_gyre']      = pd.Series(fh.variables['BSF_gyre'],index=time)
    dfs['KE_total_SO']   = pd.Series(fh.variables['KE_total_SO'],index=time)
    dfs['KE_top_SO']     = pd.Series(fh.variables['KE_top_SO'],index=time)
    dfs['KE_bottom_SO']  = pd.Series(fh.variables['KE_bottom_SO'],index=time)
    dfs['OHC_total_SO']  = pd.Series(fh.variables['OHC_total_SO'],index=time)
    dfs['OHC_top_SO']    = pd.Series(fh.variables['OHC_top_SO'],index=time)
    dfs['OHC_bottom_SO'] = pd.Series(fh.variables['OHC_bottom_SO'],index=time)
    dfs['SOM_index']     = pd.Series(fh.variables['SOM_index'],index=time)
    
    fh.close()
    
    return dfs
