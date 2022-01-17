% This code solves the idealized model of sea ice and climate described by
% Roach et al. (2022; see reference below) in the "NoIce_NoDiff"
% configuration at a single latitude, which is described by eq. (7) in
% Roach et al. (2022):
%   cw*dT/dt = a*S - [A+B*(T-T_f)] + Fb .
% The code solves this ODE for T(t) by using an n-term Fourier decomposition
% of insolation S(t), analytically solving each term in the decomposition and
% omitting the transient spinup term, and then summing up the solutions.
%
% For code to solve the model in all configurations using numerical time
% stepping, which includes more comments about the model, see the file
%   sea_ice_EBM_R22.m at https://eisenman-group.github.io
%
% This code gives a solution for each latitude that is approximately
% equivalent to running sea_ice_EBM_R22('Config=''NoIce_NoDiff'''),
% although this code runs considerably faster by avoiding numerical time
% stepping and spinup issues.
%
% Inputs: Any changes to parameter values (see Example below).
%
% Outputs: Seasonal cycle (neglecting initial transients) of
% t (time in yrs), and T (surface temperature in deg C).
%
% Example:
% [t,~,~,T] = sea_ice_EBM_R22('Config=''NoIce_NoDiff''','A=138.6');
% lati=362; dx=1/400; lat=-asind(-dx/2+lati*dx); % lati the index of the latitude to be plotted
% [t2,T2] = sea_ice_EBM_R22_eq7(['lat=' num2str(lat)]);
% plot(t,T(:,lati),t2,T2,'--'), ylabel(['T (^oC) at ' num2str(-lat,3) '^oS']), xlabel('t (yrs)')
%
% Code by Ian Eisenman (eisenman@ucsd.edu), 2022.
%
% Reference:
% L. Roach, I. Eisenman, T.J.W. Wagner, E. Blanchard-Wrigglesworth, C. Bitz
%   (2022). Asymmetry in the seasonal cycle of Antarctic sea ice due to
%   insolation. In press, Nature Geoscience.
%
%--------------------------------------------------------------------------
function [t_c, T_c] = sea_ice_EBM_R22_eq7(varargin)
%%Model parameters --------------------------------------------------------
lat=-64.5;    %latitude to consider
Sol='SHSol';  %solar forcing ('SHSol', 'IdealSol', or 'NHSol')
A  = 138.6;   %OLR when T = T_f (W m^-2)
B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)
cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
a0 = 0.7;     %ice-free co-albedo at equator
a2 = 0.1;     %ice-free co-albedo spatial dependence
Fb = 5;       %heat flux from ocean below (W m^-2)
s0 = 420;     %IdealSol insolation at equator  (W m^-2)
s1 = 338;     %IdealSol insolation seasonal dependence (W m^-2)
s2 = 240;     %IdealSol insolation spatial dependence (W m^-2)
trunc=50;     %number of Fourier coefficients to retain in solution (max is nt/2)
nt = 1e3;     %# of timesteps per year (for insolation calculation)
%--------------------------------------------------------------------------
%parameter value changes as input to function (batch mode)
if nargin>0, for j=1:length(varargin), eval([varargin{j} ';']), end, end
%--------------------------------------------------------------------------
dt = 1/nt;
t = dt:dt:1;       %time array (yrs)
x=-sind(lat);
SHSol=daily_insolation_subroutine(lat,t*365.24);
NHSol=daily_insolation_subroutine(-lat,t*365.24);
IdealSol=s0-s1*cos(2*pi*t)*(-x)-s2*x^2;
if strcmp(Sol,'SHSol'), S=SHSol; end
if strcmp(Sol,'IdealSol'), S=IdealSol; end
if strcmp(Sol,'NHSol'), S=NHSol; end

% Compute Fourier Series decomposition of insolation (using fft)
Sf=fft(S)/nt;
ak0=2*Sf(1); % coefficient a_0 (twice the time average)
ak=[real(Sf(end:-1:nt/2+2)+Sf(2:nt/2)) Sf(nt/2+1)]; % cosine coefficients a_k
bk=[imag(Sf(end:-1:nt/2+2)-Sf(2:nt/2)) 0]; % sine coefficients b_k
% Analytical solution for each Fourier Component
% dx/dt = c*cos(k*t)-a-b*x ==> x(t) = -a/b+c/(b^2+k^2)*(b*cos(k*t)+k*sin(k*t))
% dx/dt = c*sin(k*t)-a-b*x ==> x(t) = -a/b+c/(b^2+k^2)*(b*sin(k*t)-k*cos(k*t))
a=a0-a2*sin(deg2rad(lat))^2; % coalbedo
T=zeros(1,nt)+(Fb-A+a*ak0/2)/B;
for k=1:trunc % sum up analytical solutions
    T = T + (a*cw)/(B^2+(2*pi*k*cw)^2) * ( (ak(k)*B/cw-bk(k)*2*pi*k)*cos(2*pi*k*t) ...
        + (ak(k)*2*pi*k+bk(k)*B/cw)*sin(2*pi*k*t));
end

if nargout>0 %save output
    t_c=t';
    T_c=T';
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
