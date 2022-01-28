% This code numerically solves the idealized model of sea ice and climate
% described by Roach et al. (2022; see reference below).
%
% This is based on the model of Wagner & Eisenman (2015, hereafter WE15),
% with the main difference being that this model uses realistic insolation
% whereas WE15 used an idealized sinusoidal insolation field. More comments
% regarding the equations in this code can be found in the WE15 model code:
%   sea_ice_EBM_WE15.m at https://eisenman-group.github.io
%
% The insolation field is calculated using the code of Huybers and Eisenman
% (2006), which is included in this file as a subroutine. A more complete
% version of this code (including past orbital parameter values and more
% comments) is available is a separate file:
%   daily_insolation.m at https://eisenman-group.github.io
%
% The results presented in Roach et al. (2022) were generated using python
% code that is available at the link given at the end of the paper. This
% Matlab code produces results that have slight numerical differences
% compared with the python results.
%
% Inputs: Any changes to parameter values (see Example below).
% 
% Outputs: Climatological mean seasonal cycle during final 20 years of 
% t (time in yrs), iceA (ice area in 10^6 km^2), E (enthalpy in J/m^2), 
% and T (surface temperature in deg C). If no outputs are specified in the
% function call, then the fields are saved to sea_ice_EBM_R22.mat .
%
% Example:
% [t1,iceA1] = sea_ice_EBM_R22;
% [t2,iceA2] = sea_ice_EBM_R22('Sol=''IdealSol''','A=195.0');
% plot(t1,iceA1,t2,iceA2), ylabel('Ice area (10^6 km^2)'), xlabel('Time (yrs)')
%
% This code was compiled for sharing based on the model described by Roach
% et al. (2022), drawing on code from the github.io links above, by 
% Ian Eisenman (eisenman@ucsd.edu), 2022.
%
% References:
% L. Roach, I. Eisenman, T.J.W. Wagner, E. Blanchard-Wrigglesworth, C. Bitz
%   (2022). Asymmetry in the seasonal cycle of Antarctic sea ice due to
%   insolation. In press, Nature Geoscience.
% T.J.W. Wagner and I. Eisenman (2015). How climate model complexity
%   influences sea ice stability. J Climate 28, 3998-4014.
% P. Huybers and I. Eisenman (2006). Integrated summer insolation
%   calculations. NOAA/NCDC Paleoclimatology Program, Data Contribution
%   Series #2006-079.
%
%--------------------------------------------------------------------------
function [t_c, iceA_c, E_c, T_c] = sea_ice_EBM_R22(varargin)
%%Model parameters --------------------------------------------------------
Config='Full';%model configuration ('Full', 'NoIce', or 'NoIce_NoDiff')
Sol='SHSol';  %solar forcing ('SHSol', 'IdealSol', or 'NHSol')
D  = 0.625;   %diffusivity for heat transport (W m^-2 K^-1)
A  = 194.8;   %OLR when T = T_f (W m^-2)
B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)
cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
a0 = 0.7;     %ice-free co-albedo at equator
a2 = 0.1;     %ice-free co-albedo spatial dependence
ai = 0.4;     %co-albedo where there is sea ice
Fb = 5;       %heat flux from ocean below (W m^-2)
k  = 2;       %sea ice thermal conductivity (W m^-2 K^-1)
Lf = 9.5;     %sea ice latent heat of fusion (W yr m^-3)
s0 = 420;     %IdealSol insolation at equator  (W m^-2)
s1 = 338;     %IdealSol insolation seasonal dependence (W m^-2)
s2 = 240;     %IdealSol insolation spatial dependence (W m^-2)
cg = 0.098;   %ghost layer heat capacity(W yr m^-2 K^-1)
tau = 3e-5;   %ghost layer coupling timescale (yr) % XX
%%time-stepping parameters ------------------------------------------------
n  = 400;     %# of evenly spaced latitudinal gridboxes (equator to pole)
nt = 1e3;     %# of timesteps per year (limited by numerical stability)
dur= 100;     %# of years for the whole run
%--------------------------------------------------------------------------
%parameter value changes as input to function (batch mode)
if nargin>0, for j=1:length(varargin), eval([varargin{j} ';']), end, end
%--------------------------------------------------------------------------
%%NoDiff switch -----------------------------------------------------------
if strcmp(Config,'NoIce_NoDiff'), D=0; end
%%Grid --------------------------------------------------------------------
dx = 1/n;               %grid box width
x = (dx/2:dx:1-dx/2)';  %grid
dt = 1/nt;
%%Diffusion Operator ------------------------------------------------------
xb = (dx:dx:1.0-dx)';
lambda=D/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;
diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);
%%Definitions for implicit scheme for Tg
cg_tau = cg/tau;
dt_tau = dt/tau;
dc = dt_tau*cg_tau;
kappa = (1+dt_tau)*eye(n)-dt*diffop/cg;
%%Insolation
ty = dt:dt:1;
xx=repmat(-x,[1,length(ty)]); tt=repmat(ty,[length(x),1]);
lat=rad2deg(asin(xx)); day=tt*365.24;
SHSol=daily_insolation_subroutine(lat,day);
NHSol=daily_insolation_subroutine(-lat,day);
IdealSol=s0-s1*cos(2*pi*tt).*xx-s2*xx.^2;
if strcmp(Sol,'SHSol'), S=SHSol; end
if strcmp(Sol,'IdealSol'), S=IdealSol; end
if strcmp(Sol,'NHSol'), S=NHSol; end
S=[S S(:,1)];
%%Further definitions -----------------------------------------------------
M = B+cg_tau;
aw= a0-a2*x.^2;   % ice-free albedo
kLf = k*Lf;
%%Initial condition -------------------------------------------------------
T = 7.5+20*(1-2*x.^2);
Tg = T; E = cw*T;
%%Set up output arrays, saving every timestep
E1000 = nan(n,dur*nt); T1000 = E1000;
%%Integration (see WE15_NumericIntegration.pdf)----------------------------
% Loop over Years ---------------------------------------------------------
if strcmp(Config,'NoIce') || strcmp(Config,'NoIce_NoDiff')
    ice=0;
else
    ice=1;
end
disp('starting run')
for years = 1:dur
    % Loop within One Year-------------------------------------------------
    for i = 1:nt
        %if mod(i,nt/100)==0 %store 100 timesteps per year
        E1000(:,(years-1)*nt+i) = E;
        T1000(:,(years-1)*nt+i) = T;
        %end
        % forcing
        if ice
            alpha = aw.*(E>0) + ai*(E<0);        % WE15 Eq. (4)
        else
            alpha = aw;
        end
        C = alpha.*S(:,i) + cg_tau*Tg-A;
        % surface temperature
        if ice
            T0 =  C./(M-kLf./E);                 %WE15 Eq. (A3)
            T = E/cw.*(E>=0)+T0.*(E<0).*(T0<0);  %WE15 Eq. (9)
            % Forward Euler for E
            E = E+dt*(C-M*T+Fb);                 %WE15 Eq. (A2)
            % Implicit Euler for Tg
            Tg = (kappa-diag(dc./(M-kLf./E).*(T0<0).*(E<0)))\ ...
                (Tg + (dt_tau*(E/cw.*(E>=0)+(ai*S(:,i+1) ...
                -A)./(M-kLf./E).*(T0<0).*(E<0))));        %WE15 Eq. (A1)
        else
            T = E/cw;
            E = E+dt*(C-M*T+Fb);
            Tg = kappa\ (Tg + (dt_tau*(E/cw)));
        end
    end
    if mod(years,25)==0, disp(['year ' num2str(years) ' complete']), end
end
% -------------------------------------------------------------------------
%output mean climatology during final 20 years
Ec=squeeze(mean(reshape(E1000(:,end-20*nt+1:end)',nt,[],n),2)); % climatological mean
Tc=squeeze(mean(reshape(T1000(:,end-20*nt+1:end)',nt,[],n),2)); % climatological mean
tc = ty(1:nt)';
% ice area: interpolate to E=0
iceA=nan(size(tc));
earth_A=510; % surface area of earth in 10^6 km^2
for m=1:length(iceA)
    iceA(m)=earth_A/2*(1-interp1(Ec(m,:),x,0));
    if min(Ec(m,:))>0, iceA(m)=0; end
end

if nargout>0 %save output
    t_c=tc;
    iceA_c=iceA;
    E_c=Ec;
    T_c=Tc;
else
    save sea_ice_EBM_R22.mat tc iceA Ec Tc
    disp('output saved to sea_ice_EBM_R22.mat')
end


function Fsw = daily_insolation_subroutine(lat,day)
%
% Description:
%   Shortened version of daily_insolation.m at https://eisenman-group.github.io
%   Computes daily average insolation in W/m2 as a function of day and
%   latitude using modern orbital parameters.
%
% Detailed description of calculation:
%   If using calendar days, solar longitude is found using an approximate
%   solution to the differential equation representing conservation of
%   angular momentum (Kepler's Second Law).  Given the orbital parameters
%   and solar longitude, daily average insolation is calculated exactly
%   following Berger 1978.
%
% References:
%   Berger A. (1978). Long-term variations of daily insolation and
%     Quaternary climatic changes. Journal of Atmospheric Science, 35(12),
%     2362-2367.
%
% Authors:
%   Ian Eisenman and Peter Huybers, Harvard University, 2006
%   eisenman@post.harvard.edu
%

% === Get orbital parameters ===
m=[0 0.017236 101.37 23.446];
ecc=m(:,2); % eccentricity
% add 180 degrees to omega (see lambda definition, Berger 1978 Appendix)
omega=m(:,3)+180; % longitude of perihelion (precession angle)
omega=unwrap(omega*pi/180); % remove discontinuities (360 degree jumps)
epsilon=m(:,4)*pi/180; % obliquity angle

% === Calculate insolation ===
lat=lat*pi/180; % latitude
% lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)

% use calendar days
% estimate lambda from calendar day using an approximation from Berger 1978 section 3
delta_lambda_m=(day-80)*2*pi/365.2422;
beta=(1-ecc.^2).^(1/2);
lambda_m0=-2*( (1/2*ecc+1/8*ecc.^3).*(1+beta).*sin(-omega)-...
    1/4*ecc.^2.*(1/2+beta).*sin(-2*omega)+1/8*ecc.^3.*(1/3+beta).*(sin(-3*omega)) );
lambda_m=lambda_m0+delta_lambda_m;
lambda=lambda_m+(2*ecc-1/4*ecc.^3).*sin(lambda_m-omega)+...
    (5/4)*ecc.^2.*sin(2*(lambda_m-omega))+(13/12)*ecc.^3.*sin(3*(lambda_m-omega));
So=1365; % solar constant (W/m^2)
delta=asin(sin(epsilon).*sin(lambda)); % declination of the sun
Ho=acos(-tan(lat).*tan(delta)); % hour angle at sunrise/sunset
% no sunrise or no sunset: Berger 1978 eqn (8),(9)
Ho( ( abs(lat) >= pi/2 - abs(delta) ) & ( lat.*delta > 0 ) )=pi;
Ho( ( abs(lat) >= pi/2 - abs(delta) ) & ( lat.*delta <= 0 ) )=0;

% Insolation: Berger 1978 eq (10)
Fsw=So/pi*(1+ecc.*cos(lambda-omega)).^2 ./ ...
    (1-ecc.^2).^2 .* ...
    ( Ho.*sin(lat).*sin(delta) + cos(lat).*cos(delta).*sin(Ho) );
