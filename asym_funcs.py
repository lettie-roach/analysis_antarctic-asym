import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
import os
xr.set_options(keep_attrs=True)
from scipy.fft import fft, ifft
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':14})
rc('text', usetex=True)

a_earth = 510072000./1e6 # million km^2

firstofmonthind = [1,32,60,91,121,152,182,213,244,274,305,335]
firstofmonthlabel = ['01-Jan','01-Feb','01-Mar','01-Apr','01-May','01-Jun','01-Jul',
                    '01-Aug','01-Sep','01-Oct','01-Nov','01-Dec']

shfirstofmonthind = [1,32,63,93,124,154,185,216,244,275,305,335]
shfirstofmonthlabel = ['01-Jul','01-Aug','01-Sep','01-Oct','01-Nov','01-Dec','01-Jan','01-Feb','01-Mar','01-Apr','01-May','01-Jun']


def get_annual_harmonic (xdata, ydata):
    """
    Compute the annual harmonic using FFT
    """
    
    N = len(ydata)
    yinv = fft(ydata)/N
    
    annual_harmonic = yinv[0] + np.real(yinv[N-1] + yinv[1])*np.cos(2.*np.pi*xdata) + np.imag(yinv[N-1] - yinv[1])*np.sin(2.*np.pi*xdata)
    
    return annual_harmonic
                                                                                                                          
                                                                                                                          
def _nanlinregress(x, y):
    """
    Calls scipy linregress only on finite numbers of x and y
    """
    
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        # empty arrays passed to linreg raise ValueError:
        # force returning an object with nans:
        return scipy.stats.linregress([np.nan], [np.nan])
    return scipy.stats.linregress(x[finite], y[finite])



def linregress(first_samples, second_samples, dim):
    """
    Apply the _nanlinregress function to an xarray Dataset
    """
    slope, intercept, r_value, p_value, std_err = xr.apply_ufunc(_nanlinregress,
                       first_samples, second_samples,
                       input_core_dims  = [[dim], [dim]], 
                       output_core_dims = [[],[],[],[],[]],
                       vectorize=True)    
    return slope, intercept, r_value, p_value, std_err






def asym_1d (mydata):
    """
    Compute the difference between the length of the growth and retreat seasons
    Seasons are defined time of maximum and minimum occurrence
    
    Note that this only works for extrema that are within the time axis - 
     i.e. not if you have an extrema near the end and beginning of the period
    """
    dmin =  np.argwhere(mydata<=np.nanmin(mydata)).mean()
    dmax =  np.argwhere(mydata>=np.nanmax(mydata)).mean()
    nt = len(mydata)

    if dmax>dmin:
        growth = dmax - dmin
        retreat = nt - growth
    else:
        retreat = dmin - dmax
        growth = nt - retreat
        
    asym = growth - retreat

    return asym, dmin, dmax


def asym_xr (ds, dim) :                                                                                                                     
    """
    Apply the asym_1d function to an xarray Dataset
    """ 
    asym, dmin, dmax = xr.apply_ufunc(asym_1d,
                   ds, 
                   input_core_dims  = [[dim]], 
                   output_core_dims = [[],[],[]],
                   vectorize=True)

    
    return asym, dmin, dmax
                                                                                                                          
                                                                                                                          
def switch_nh(siads):
    """
    Shift NH data by half a year in a dataset that contains NH and SH data
    """ 
    nhdays = np.append(np.arange(183,366,1),np.arange(1,183,1))
    nhds = siads[[f for f in siads.variables if 'nh' in f]]
    nhds = nhds.sel(day=nhdays)
    nhds['day'] = np.arange(1,366,1)
    shds = siads[[f for f in siads.variables if 'sh' in f]]
    siads = xr.merge([nhds,shds])
    
    return siads


def switch_sh(siads):
    """
    Shift SH data by half a year in a dataset that contains NH and SH data
    """ 
    shdays = np.append(np.arange(182,366,1),np.arange(1,182,1))
    shds = siads[[f for f in siads.variables if 'sh' in f]]
    shds = shds.sel(day=shdays)
    shds['day'] = np.arange(1,366,1)
    nhds = siads[[f for f in siads.variables if 'nh' in f]]
    siads = xr.merge([nhds,shds])
    
    return siads


def grid_area_regll(lat,lon):
    """
    Compute the area of grid cells on a regular lat lon grid, given lats and lons
    """ 
    to_rad = 2. *np.pi/360.
    r_earth = 6371.22 # km
    con = r_earth*to_rad
    clat = np.cos(lat*to_rad)
    dlon = lon[2] - lon[1]
    dlat = lat[2] - lat[1]
    dx = con*dlon*clat
    dy = con*dlat
    dxdy = dy*dx
    garea = np.swapaxes(np.tile(dxdy,(len(lon),1)),0,1)
    latl = np.swapaxes(np.tile(lat,(len(lon),1)),0,1)
    nh_area = np.where(latl<0.,0.,garea)
    sh_area = np.where(latl>0.,0.,garea)
    
    return garea, nh_area, sh_area




def we15_model (exptype = 'Full',nspace=400, nts=1000, mysolar=0., mycw=9.8, myA=193, dur=100, myD=0.6, myFb=4):
    """
    Run the WE15 model. Can vary parameters, resolution, run length and complexity of model as input to this function
    """                                                                                                                        
                                                                                                                          

    # Physical parameters
    D  = myD #0.6  # diffusivity for heat transport (W m^-2 K^-1)
    S1 = 338.;     # insolation seasonal dependence (W m^-2)
    A  = myA #193  # OLR when T = 0 (W m^-2)
    B  = 2.1       # OLR temperature dependence (W m^-2 K^-1)
    cw = mycw      # ocean mixed layer heat capacity (W yr m^-2 K^-1)
    S0 = 420.      # insolation at equator  (W m^-2)
    S2 = 240.      # insolation spatial dependence (W m^-2)
    a0 = 0.7       # ice-free co-albedo at equator
    a2 = 0.1       # ice=free co-albedo spatial dependence
    ai = 0.4       # co-albedo where there is sea ice
    Fb = myFb #4   # heat flux from ocean below (W m^-2)
    k  = 2;        # sea ice thermal conductivity (W m^-2 K^-1)
    Lf = 9.5;      # sea ice latent heat of fusion (W yr m^-3)
    cg = 0.01*cw;  # ghost layer heat capacity(W yr m^-2 K^-1)
    tau = 3e-5;    # ghost layer coupling timescale (yr)
    a_earth = 510072000./1e6 # million km^2

    # Time stepping parameters
    ##The default run in WE15, Fig 2 uses the time-stepping parameters: -------
    #n=400; % # of evenly spaced latitudinal gridboxes (equator to pole)
    #nt=1e3; % # of timesteps per year (approx lower limit of stability) 
    #dur=200; % # of years for the whole run

    n  = nspace;
    nt = nts;
    dt = 1/nt;

    
    #Spatial Grid -------------------------------------------------------------
    dx = 1.0/n    #grid box width
    x = np.arange(dx/2,1+dx/2,dx) #native grid
    xb = np.arange(dx,1,dx)

    ##Diffusion Operator (WE15, Appendix A) -----------------------------------
    lam = D/dx**2*(1-xb**2)
    L1=np.append(0, -lam); L2=np.append(-lam, 0); L3=-L1-L2
    diffop = - np.diag(L3) - np.diag(L2[:n-1],1) - np.diag(L1[1:n],-1);

    ##Definitions for implicit scheme on Tg
    cg_tau = cg/tau; 
    dt_tau = dt/tau; 
    dc = dt_tau*cg_tau;
    kappa = (1+dt_tau)*np.identity(n)-dt*diffop/cg;

    ##Seasonal forcing (WE15 eq.3)
    ty = np.arange(dt/2,1+dt/2,dt)
    if np.all(mysolar==0.): # if an insolation field is not provided, use the idealized function
        S = (np.tile(S0-S2*x**2,[nt,1])- np.tile(S1*np.cos(2*np.pi*ty),[n,1]).T*np.tile(x,[nt,1])); #totally symmetric at all lats
    else:
        S = mysolar
    
    ##Further definitions
    M = B+cg_tau; 
    aw = a0-a2*x**2      # open water albedo
    kLf = k*Lf;
    
    #Set up output arrays, saving all timesteps of all years
    Efin  = np.zeros([dur,nt,n])
    Tfin  = np.zeros([dur,nt,n])
   
    #Initial conditions ------------------------------------------------------
    T = 7.5+20*(1-2*x**2);
    Tg = T; E = cw*T;
    
    #Integration (see WE15_NumericIntegration.pdf)----------------------------
    #Loop over Years ---------------------------------------------------------
    sia = np.zeros([dur,nt])
    for years in range(0,dur):
        #Loop within One Year-------------------------------------------------
        for i in range(0,int(nt)):
            #store spatial fields for all years

            #forcing
            alpha = aw*(E>0) + ai*(E<0)     #WE15, eq.4
            
            if 'NoAlb' in exptype:
                alpha = aw

            C = alpha*S[i,:]+cg_tau*Tg-A
            
            #surface temperature
            T0 = C/(M-kLf/E)                 #WE15, eq.A3
                           
            if 'NoTherm' in exptype:
                T = E/cw
                E = E+dt*(C-M*T+Fb);                 #WE15, eq.A2, Forward Euler on E
                Tg = np.linalg.solve(kappa,Tg+(dt_tau*(E/cw)))
                
                       
            else:
                T = E/cw*(E>=0)+T0*(E<0)*(T0<0);  #WE15, eq.9
                E = E+dt*(C-M*T+Fb);                 #WE15, eq.A2, Forward Euler on E
                Tg = np.linalg.solve(kappa-np.diag(dc/(M-kLf/E)*(T0<0)*(E<0)),
                                     Tg+(dt_tau*(E/cw*(E>=0)+(ai*S[i,:]-A)/(M-kLf/E)*(T0<0)*(E<0)))) #Implicit Euler on Tg
                
            Efin[years,i,:] = E
            Tfin[years,i,:] = T
                   
    lat = np.arcsin(np.linspace(0,1,nspace))*180./np.pi
    E_all = xr.DataArray(Efin,dims=('year','day','lat'),coords = {'year':np.arange(1,dur+1,1), 'day':np.linspace(1,366,nts), 'lat':lat}).to_dataset(name='E')
    T_all = xr.DataArray(Tfin,dims=('year','day','lat'),coords = {'year':np.arange(1,dur+1,1), 'day':np.linspace(1,366,nts), 'lat':lat}).to_dataset(name='T')
    
    ds = xr.merge([E_all, T_all]) 
    
    return ds







