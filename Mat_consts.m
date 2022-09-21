function [mu0, Msat, rhok1, D, Halpha, Htheta] = Mat_consts()


% Universal Constant (permeability of free space)
mu0 = 4*pi*1e-7;

% Material Constants
% Come from manufacturer
Msat = 544296;              % magnetic saturation (Amp/m)
% rhok1 = 0.2e5;              % anisotropy constant calibrated (Joule/m^3)
% rhok1 = 1.9e5;              % anisotropy constant given by manufacturer
% rhok1 = 0.05e5;              % anisotropy constant exact match
rhok1 = 0.89e5;
Halpha = 23000; % magnetic domain saturation  (Amp/m)
% Htheta=510000;  % magnetic rotation saturation  (Amp/m)
Htheta=(2*rhok1)/(mu0*Msat);
% Htheta=555569.7486;

% Demagnetization Factors
% This is the demagnetization from the closed form solution
% http://www.magpar.net/static/magpar/doc/html/demagcalc.html
D=[ 0.06741246546340422,0,0;
    0,0.4662937598818836,0;
    0,0,0.4662937598818836];