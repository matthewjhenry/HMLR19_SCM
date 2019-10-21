% Taken from the code accompanying "Analytic radiative-advective equilibrium as a model for high-latitude climate" by Cronin & Jansen (2016, GRL)

function [ T, TS, olr , tau, I, F] = grae_w( p, tau0, D, FS, FA, alpha, R_cp, b, beta)
%grae: gray radiative-advective equilibrium with a window 
% function calculates temperature profile as a function of pressure for
% 'radiative-advective equilibrium'
%
% assumes power-law dependence of optical thickness on p, \tau = \tau_0
% (p/p0)^2
%
% and assumes that atmospheric heating FA is distributed in \tau according
% to Q_net = FA*(1-(\tau/\tau_0)^b) (if b=1 this is uniform, b=1/2 is
% uniform in pressure)
% surface heating FS is equal to absorbed SW at surface, plus subsurface
% heat flux; atmospheric heating FA is equal to the atmospheric absorbed
% solar radiation, plus the moist static energy flux convergence in the
% column from large-scale dynamics
%
% input arguments:
%    p: array of pressures, also output; units = hPa, max(p) is assumed to
%      be surface
%    tau0: total gray IR optical thickness of atmosphere (unitless)
%    D: radiation diffusivity factor (unitless)
%    FS: solar heat flux absorbed at surface, plus subsurface upward heat
%      flux (units W/m^2)
%    FA: solar heat flux absorbed by atmosphere, plus large-scale moist
%      static energy flux convergence (units W/m^2)
%    alpha: parameter in neutral convective lapse rate, 
%      d(log T)/d(log p) = alpha*R_cp (unitless)
%    R_cp: ratio of molar gas constant to specific heat at constant
%      pressure (unitless)
%    b: power for distribution of heating
%    beta: spectral width of atmospheric window (0-1)
%    lin : 0 for nonlinear LW radiation and 1 for linearized LW radiation
% 
%
% output arguments:
%    p: output array of pressures (hPa)
%    T: array of temperatures in radiative-advective equilibrium (K) --
%      T(end) is surface temperature
%    stability: 1 if atmosphere is convectively stable everywhere; 
%              -1 if unstable at surface but not free troposphere
%              -2 if unstable in free troposphere but not at surface
%              -3 if unstable both in free troposphere and surface
%      (based on comparing to dry-adiabatic d(ln T)/d(ln p) = R/cp = 2/7 for a diatomic gas

% constants
sigb = 5.67e-8;

% option to make a figure or not
do_figure = 1;

%normalized pressure
pn   = p./max(p);

% midpoint pressures
pmid = (p(1:end-1) + p(2:end))/2;

tau=tau0*pn.^2;

% surface temperature generally differs from T(tau0)
TS = ((1/sigb)*(FS*(1 + D*tau0/2)+FA/2*(1 + D*tau0*(1-1/(b+1))))/(1+beta*D*tau0/2)).^(1/4);

% general temperature expressions
T4 = (1/(2*sigb*(1-beta)))*((FS-beta*sigb*TS.^4)*(1+D*tau0.*pn.^2)+FA*(b/(D*tau0).*pn.^(2*(b-1)) + 1 + D*tau0.*(pn.^2).*(1-pn.^(2*b)/(b+1))));

T = T4.^(1/4);

olr = integ_fr(TS,T,tau0,p);

I = (FS+FA-beta*sigb*TS.^4).*(1+D*tau0.*pn.^2) - FA * tau0./(1+b).*pn.^(2*(b+1));

F = FS + FA.*(1-pn.^(2*b));

end

