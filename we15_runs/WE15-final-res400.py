import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, '/glade/u/home/lettier/analysis/analysis_antarctic-asym/pub/')
import asym_funcs as af
import time

start = time.time()

mydir = './we15_final_runs/'

# regular parameter settings
D = 0.625
Fb = 5.
nspace = 400

ds = xr.open_dataset('daily_insolation_n'+str(nspace)+'_nh.nc')
nh = ds.S.values


A = 194.8
ds = af.we15_model(exptype='WE15',nspace=nspace, dur=100, mysolar=nh, myD=D, myA=A, myFb=Fb)
myname = 'WE15_SolNH_A'+str(A)+'_res'+str(nspace)
ds['names'] = myname
ds = ds.set_coords('names')
ds = ds.isel(year=slice(-20,None)).mean(dim='year')
ds.to_netcdf(mydir+myname+'.nc')	 



ds = xr.open_dataset('daily_insolation_n'+str(nspace)+'_sh.nc')
sh = ds.S.values

for s, sol in enumerate(['SH','Ideal']):
	if sol=='SH':
		mysol = sh
	else:	
		mysol = 0.

	A = [194.8,195.][s]
	ds = af.we15_model(exptype='WE15',nspace=nspace, dur=100, mysolar=mysol, myD=D, myA=A, myFb=Fb)
	myname = 'WE15_Sol'+sol+'_A'+str(A)+'_res'+str(nspace)
	ds['names'] = myname
	ds = ds.set_coords('names')
	ds = ds.isel(year=slice(-20,None)).mean(dim='year')
	ds.to_netcdf(mydir+myname+'.nc')	 

	A = [201.6, 201.5][s]
	ds = af.we15_model(exptype='WE15-NoAlbNoTherm',nspace=nspace, dur=100, mysolar=mysol, myD=D, myA=A, myFb=Fb)
	myname = 'WE15_NoAlbNoTherm_Sol'+sol+'_A'+str(A)+'_res'+str(nspace)
	ds['names'] = myname
	ds = ds.set_coords('names')
	ds = ds.isel(year=slice(-20,None)).mean(dim='year')
	ds.to_netcdf(mydir+myname+'.nc')	 

	A = [138.6, 143.6][s]
	ds = af.we15_model(exptype='WE15-NoAlbNoTherm',nspace=nspace, dur=100, mysolar=mysol, myD=0., myA=A, myFb=Fb)
	myname = 'WE15_NoAlbNoThermNoDiff_Sol'+sol+'_A'+str(A)+'_res'+str(nspace)
	ds['names'] = myname
	ds = ds.set_coords('names')
	ds = ds.isel(year=slice(-20,None)).mean(dim='year')
	ds.to_netcdf(mydir+myname+'.nc')	 
