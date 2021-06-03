function Bc = morris_biocoeffs(sm)
%
%-------function help------------------------------------------------------
% NAME
%   morris_biocoeffs.m
% PURPOSE
%   Calculate the biomass coefficients for Bc=aD+bD^2+c where D is
%   the depth at peak biomass
% USAGE
%   Bc = morris_biocoeffs(sm)
% INPUTS
%   sm - Saltmarsh instance or struct of Saltmarsh properties
%        NumSpecies - number of species
%        MinSpDepth - minimumdepth for each species
%        MaxSpDepth - maximum depth for each species (m)
%        MaxBiomass - maximum biomass for each species (kg/m2)
%        SettlingAlpha - coefficient for biomass dependent enhanced settling rate (m/s)
%        SettlingBeta - exponent for biomass dependent enhanced settling rate  
% OUTPUTS
%   Bc - matrix (n x 3) of biomass coefficients for each species
% SEE AlSLO
%   see Saltmarsh class in Asmita
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
% 
    nsp = sm.NumSpecies;
    dmn = sm.MinSpDepth;
    dmx = sm.MaxSpDepth;
    bmx = sm.MaxBiomass;
    Bc = zeros(nsp,3);
    dpk = (dmx+dmn)/2;            %peak biomass depth
    for k = 1:nsp
        vb = [bmx(k); 0; 0];
        mD = [dpk(k) dpk(k)^2 1; dmx(k) dmx(k)^2 1; dmn(k) dmn(k)^2 1];
        Bc(k,:) = (mD)\vb;
    end         
end