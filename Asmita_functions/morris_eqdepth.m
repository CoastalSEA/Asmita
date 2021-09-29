function deq = morris_eqdepth(sm,cn,qm,dslr)
%
%-------function help------------------------------------------------------
% NAME
%   morris_eqdepth.m
% PURPOSE
%   Solve the Morris equation for equilibrium depth
% USAGE
%   deq = morris_eqdepth(sm,cn,qm,dslr)
% INPUTS
%   sm - Saltmarsh instance or struct of Saltmarsh properties, with:
%        NumSpecies - number of species
%        MinSpDepth - minimum depth for each species (m)
%        MaxSpDepth - maximum depth for each species (m)
%        MaxBiomass - maximum biomass for each species (kg/m2)
%        SpeciesProduct - species productivity (m2/kg/yr)
%        SettlingAlpha - coefficient for biomass dependent enhanced settling rate (m/s)
%        SettlingBeta - exponent for biomass dependent enhanced settling rate  
%   cn - struct of abbreviated Constants values (eg cn.g for Gravity)
%   qm - sediment loading (s^-1)
%   dslr - rate of sea level rise (m/s)
% OUTPUTS
%   deq - equilbirum depth of marsh
% SEE AlSLO
%   see Saltmarsh class in Asmita
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
%
    kbm = sm.SpeciesProduct';  
    kbm = kbm/cn.y2s;           %rate of biomass production (m^2/kg/s)
    %kbm = kbm/acb;             %correct for bed density
    Bc = morris_biocoeffs(sm);  %biomass coefficients
    %
    KBM = repmat(kbm,1,3);      %replicate kbm for each coefficient
    kbc = Bc.*KBM;              %multiply each coefficient by kbm
    cfs = sum(kbc,1);           %sum each coefficient for all species
    cfs(3)=cfs(3)+qm;           %add in the sediment loading
    coefs = [cfs -dslr];        %construct polynomial
    coefs = coefs(:,[2 1 3 4]); %swap coefs a and b to give
                                %f(D)=kbD^3+kaD^2+(qm+kc)D-slr 
                                %where kb=kbm*b summed for all species, etc.                              
    deq = roots(coefs);         %find roots
    deq = min(deq(deq>0));      %select minimum root>0
    deq = deq*isreal(deq);      %only return real value
end