function wsm = bioenhancedsettling(sm,depth,aws)
%
%-------function help------------------------------------------------------
% NAME
%   bioenhancedsettling.m
% PURPOSE
%   Adjust the vertical exchange due to biological enhanced settling
%   for the specified depth
% USAGE
%   wsm = bioenhancedsettling(sm,depth,aws)
% INPUTS
%   sm - Saltmarsh instance or struct of Saltmarsh properties
%        NumSpecies - number of species
%        MinSpDepth - minimum depth for each species (m)
%        MaxSpDepth - maximum depth for each species (m)
%        MaxBiomass - maximum biomass for each species (kg/m2)
%        SettlingAlpha - coefficient for biomass dependent enhanced settling rate (m/s)
%        SettlingBeta - exponent for biomass dependent enhanced settling rate
%   depth - depth over flat or marsh (m)
%   aws(1) = vertical exchange of element with no biology (m/s)
%   aws(2) = value for depths greater than dmx (optional, default=aws(1))   
% OUTPUTS
%   wsm - biologically enhanced settling rate
% SEE AlSLO
%   see Saltmarsh class in Asmita
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
% 
    if length(aws)<2
        aws(2) = aws(1);
    end
    dmx = sm.MaxSpDepth;       %maximum depth for each species (m)
    bmx = sm.MaxBiomass;       %maximum biomass for each species (kg/m2)
    wsa = sm.SettlingAlpha;    %coefficient for biomass dependent enhanced settling rate (m/s)
    wsb = sm.SettlingBeta;     %exponent for biomass dependent enhanced settling rate
    Bc = morris_biocoeffs(sm); %biomass coefficients
    if depth<max(dmx)
        dd  = [depth depth.^2 1];
        Bmd = Bc*dd';          %biomass for specified depth
        Bmd = Bmd.*(Bmd>0);    %only use positive values
        %------------------------------------------------------
        % Enhanced settling rate
        %------------------------------------------------------
        relBm = mean(Bmd'/(diag(bmx)),2);   %average of (Bm/bmx)
        if ~isnan(relBm) && relBm>0
            wsm = aws(1)+wsa*depth*relBm^wsb;
        else
            wsm = aws(1);
        end
    else
        wsm = aws(2);
    end
end  