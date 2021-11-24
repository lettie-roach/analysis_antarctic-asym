import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, '/glade/u/home/lettier/analysis/analysis_antarctic-asym/pub/')
import asym_funcs as af
import time

start = time.time()

mydir = './we15_cw_runs/'

# regular parameter settings
D = 0.625
Fb = 5.
nspace = 400
A = 194.8

ds = xr.open_dataset('daily_insolation_n'+str(nspace)+'_sh.nc')
sh = ds.S.values


#for cw in [6.5,9.8,13.0,16.3]:
for cw in [19.5]:
	ds = af.we15_model(exptype='WE15',nspace=nspace, dur=100, mysolar=sh, myD=D, myA=A, myFb=Fb, mycw=cw)
	myname = 'WE15_SolSH_A'+str(A)+'_D'+str(D)+'_Fb'+str(Fb)+'_Cw'+str(cw)+'_res'+str(nspace)
	ds['names'] = myname
	ds = ds.set_coords('names')
	print(ds.isel(year=slice(-20,None)))
	ds = ds.isel(year=slice(-20,None)).mean(dim='year')
	ds['Cw'] = cw
	ds.to_netcdf(mydir+myname+'.nc')	 


