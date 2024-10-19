function ws = settling_velocity(d50,g,rhow,rhos,visc,rhoc)
%
%-------function help------------------------------------------------------
% NAME
%   settling_velocity.m
% PURPOSE
%   calcualte the settling velocity using Soulsby equation
% INPUT
%   d50  = median sediment grain size (m)
%   rhoc = density of sediment(kg/m^3) - default is 0.
%   g = acceleration due to gravity (m/s2)
%   rhow = density of water (kg/m3)
%   rhos = density of sediment (kg/m3)
%   visc = kinematic viscosity of water (m2/s)
% OUTPUT
%   ws - settling velocity (m/s)
%
% Author: Ian Townend
% CoastalSEA (c)June 2015
%--------------------------------------------------------------------------
%
if nargin<6
    rhoc=0;
end
%
conc  = rhoc/rhos; %volume concentration (-)
%
% calculate settling velocity.  Mud Manual, eqn 5.7 neglecting flocculation
Ds   = d50*(g*(rhos/rhow-1)/visc^2)^(1/3);
ws   = visc/d50*(sqrt(107.3+1.049*((1-conc)^4.7)*Ds^3)-10.36);