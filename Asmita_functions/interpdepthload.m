function [depth,sedload] = interpdepthload(sm,cn,aws,qm0,dslr)
%
%-------function help------------------------------------------------------
% NAME
%   interpdepthload.m
% PURPOSE
%  Iterate to find the equilibrium depth and sediment load,qm(s^-1)
% USAGE
%   [depth,sedload] = interpdepthload(sm,cn,aws,qm0,dslr)
% INPUTS
%   sm - Saltmarsh instance or struct of Saltmarsh properties, with:
%        NumSpecies - number of species
%        MinSpDepth - minimumdepth for each species
%        MaxSpDepth - maximum depth for each species (m)
%        MaxBiomass - maximum biomass for each species (kg/m2)
%        SpeciesProduct - species productivity (m2/kg/yr)
%        SettlingAlpha - coefficient for biomass dependent enhanced settling rate (m/s)
%        SettlingBeta - exponent for biomass dependent enhanced settling rate  
%   cn - struct of abbreviated Constants values (eg cn.g for Gravity)
%   aws - vertical exchange of element with no biology (m/s)
%   qm0 - initial guess of sediment loading (s^-1)
%   dslr - rate of sea level rise (m/s)
% OUTPUTS
%   deqth - equilibrium depth (m)
%   sedload - sediment load, qm (s^-1)
% SEE ALSO
%   see Saltmarsh class in Asmita
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
%

    %concentration over marsh as a function of depth   
    cdep = sm.MarshDepthConc.Depth;
    conc = sm.MarshDepthConc.Concentration;
    %initial estimate of equilibrium depth
    depth = morris_eqdepth(sm,cn,qm0,dslr); 
    if depth<=0, depth = 0.1; end %trap no depth over marsh
    nloop = 0;
    diff = 1;
    while nloop<=10 && abs(diff)>0.001
        if depth>0
            cem = interp1q(cdep,conc,depth);
            wsm = bioenhancedsettling(sm,depth,aws);
            sedload = wsm*cem/depth;
        else
            sedload = 0;
        end
        newdepth = morris_eqdepth(sm,cn,sedload,dslr);
        diff = newdepth-depth;
        depth = depth+diff/2;
        nloop=nloop+1;
    end
    if nloop>10
        sedload= 0.0;
        depth = morris_eqdepth(sm,cn,sedload,dslr);
        sprintf('Exceed loop in Saltmarsh.EqDepthBiomass, dep=%5.3d; qm=%5.3d \n',depth,sedload)
    end
end