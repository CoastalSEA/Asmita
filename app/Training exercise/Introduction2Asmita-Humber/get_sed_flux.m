function [dvol,delV,conc] = get_sed_flux(inp,slr)
%
%-------function help------------------------------------------------------
% NAME
%   get_sed_flux.m
% PURPOSE
%   single element ASMITA model to compute amount of sedimentation
% USAGE
%   [dvol,delV] = get_sed_flux(inp,slr)
% INPUTS
%   inp - struct containing the following fields
%         inp.Volume;            %element volume (m^3)
%         inp.SurfaceArea;       %element surface area (m^2)
%         inp.Prism;             %tidal prism of channel (m^3)
%         inp.EqScaleCoeff;      %equilibrium scale coeffient, alpha
%         inp.EqShapeCoeff;      %equilibrium shape coeffient, beta 
%         inp.VerticalExchange;  %vertical exchange (m/s)
%         inp.HorizontalExchange;%horizontal exchange (m/s)
%         inp.RiverDischarge;    %river discharge (m^3/s) +ve downstream
%         inp.TransportCoeff;    %transport coefficient n (3-5)
%         inp.EqConc;            %equilibrium concentration (-)
%         inp.RiverConc;         %river load imported by advection (-)
%         inp.BedConc;           %concentration of bed (-)
%         inp.y2s;               %factor to convert from years to seconds
%   slr - rate of sea level rise (m/yr)
% OUTPUT
%   dvol - change in morphological volume in a year (m3)
%   delV - change in water volume (S*slr) in a year (m3)
%   conc - suspended sediment concentration in element (-)
% NOTE
%   if transport coefficient, n, is positive, dvol is change in water
%   volume and if n<0 it is change in sediment volume. delV is change in
%   moving surface water volume and increases with an increase is sea level
%   (ie slr>0)
% SEE ALSO
%   used in ChannelForm model (eg CF_SediData and CF_TransModel)
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    V = inp.Volume;            %element volume (m^3)  
    S = inp.SurfaceArea;       %element surface area (m^2)
    P = inp.Prism;             %tidal prism of channel (m^3)
    a = inp.EqScaleCoeff;      %equilibrium scale coeffient, alpha
    b = inp.EqShapeCoeff;      %equilibrium shape coeffient, beta 
    w = inp.VerticalExchange;  %vertical exchange (m/s)
    d = inp.HorizontalExchange;%horizontal exchange (m/s)
    q = inp.RiverDischarge;    %river discharge (m^3/s) +ve downstream
    n = inp.TransportCoeff;    %transport coefficient n (3-5)
    %
    cE = inp.EqConc;           %equilibrium concentration (-)
    cr = inp.RiverConc;        %river load imported by advection (-)
    cb = inp.BedConc;          %concentration of bed (-)
    k = cr/cE;

    C = 1/cb*(w*cE*S)/(d+q+w*S);           %constant term in dV/dt
    alpha = ((k*q+d)/(q+d))^(1/n);         %external exchange correction                    
 
    delV = S*slr;                          %water volume change (m3/yr)
    Ve = a*P^b;                            %equilbrium volume
    Gam = (alpha*Ve/(V+delV))^n;           %gamma reflects difference from equilibrium
    dvol = C*((q+d)*Gam-(k*q+d))*inp.y2s;  %morphological change (m3/yr)
    % cie = cE*Gam;                          %local equilibrium concentration in element (-)
    conc = cE*(d+q+w*S)/(w*S*Gam+k*q+d);   %sediment concentration in element (-)