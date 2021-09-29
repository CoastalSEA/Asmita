function conctable = concovermarsh(sm,wl,c0,aws,mco)
%
%-------function help------------------------------------------------------
% NAME
%   concovermarsh.m
% PURPOSE
%   Calculate time averaged concentration over the marsh for a range of
%   depths and return a dstable of depth, concentration and submergence time
%   for an emulated spring-neap cycle.
% USAGE
%   conctable = concovermarsh(sm,wl,c0,aws,mco)
% INPUT
%   sm - Saltmarsh instance or struct of Saltmarsh properties
%        NumSpecies - number of species
%        MinSpDepth - minimum depth for each species (m)
%        MaxSpDepth - maximum depth for each species (m)
%        MaxBiomass - maximum biomass for each species (kg/m2)
%        SettlingAlpha - coefficient for biomass dependent enhanced settling rate (m/s)
%        SettlingBeta - exponent for biomass dependent enhanced settling rate
%   wl - struct used in simple_tide (set unwanted constituents to zero)
%        MSL0                      mean tidel level to ordnance datum (mOD)
%        TidalAmp                  tidal elevation amplitude (m)
%        ElevPhase                 phase of elevation (ie k.x) (rads)
%        VelocityAmp               tidal velocity amplitude (m/s)
%        VelocityPhase             phase of velocity (ie k.x+phi) (rads)
%        aM2                       M2 tidal amplitufe (m)
%        aS2                       S2 tidal amplitufe (m)
%        aO1                       O1 tidal amplitufe (m)
%   c0 - background concentration adjacent to marsh (-)          
%   aws(1) = vertical exchange of element with no biology (m/s)
%   aws(2) = value for depths greater than dmx (optional, default=aws(1))   
%   mco - struct used in marsh_conc with the following settings
%        tsn  - duration of spring-neap cycle (days) [default = 14.77 days]
%        delt - time step (secs)  [default = 1 sec]
%        dmin - minimum depth used in calculations (m) [default = 0.05m];
% OUTPUTS
%   conctable - dstable of depth, concentration and submergence time
% NOTE
%   NB: result sensitive to time step and dmin
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
% 
    dm = 0;      % Initialise flags and counters
    nint = 10;
    % tidal range divided into 2*nint steps with greater resolution 
    % over the marsh if present
    dmx = sm.MaxSpDepth;                 %maximum depth for each species (m)
    if isempty(dmx)
        dmx=0;
        nno=2;
    else
        din1 = max(dmx)/nint;            %for the marsh depth range
        nno=3;
    end
    din2 = wl.TidalAmp/nint;             %for the rest of the intertidal
    dpm = zeros(nno*nint,1); cem = dpm; smt = dpm;
    tint = 1:nno*nint; 
    hw = waitbar(0,'Building Concentration Table');
    %start of depth loop-------------------------------------------
    for jd = tint       
        if dm<=max(dmx)
            dm = mco.dmin + din1*(jd-1); %depths over marsh
        else
            dm = mco.dmin + max(dmx) + din2*(jd-1-(nint*(nno==3))); %depths over flat
        end
        zm   = wl.TidalAmp-dm; %elevation of marsh(assumes z0=mtl and symmetric tide)
        %
        wsm = bioenhancedsettling(sm,dm,aws);
        %
        % Calculate the concentration using the method proposed by Krone(1987).
        [cem(jd),smt(jd)] = marsh_conc(c0,zm,wsm,wl,mco);
        dpm(jd) = dm;
        waitbar(jd/nno/nint)
    end
    close(hw)
    conctable = dstable([0;dpm], [0;cem], [0;smt],'RowNames',[0;tint'],...
                  'VariableNames',{'Depth','Concentration','Submergence'});
end