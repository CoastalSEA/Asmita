function slr = sealevelrise(time,datetime,dslr,option)
%
%-------function help------------------------------------------------------
% NAME
%   sealevelrise.m
% PURPOSE
%   Function to compute sea level rise from start of model run
% USAGE
%   slr = sealevelrise(time,datetime,dslr,option)
% INPUTS
%   time - vector of model time in seconds from t=0
%   datetime - time elapsed from Year 0 in seconds
%   dslr - rate of sea level rise (m/yr)
%   option - determine how slr is calculated
%           1 - linear slr
%           2 - exponential from 1900
%           3, etc - user can add addtional options as required
% OUTPUTS
%   slr = magnitude of sea level rise from start time to time now
%
% Author: Ian Townend
% CoastalSEA (c) Apr 2021
%--------------------------------------------------------------------------
% 
    y2s = muiConstants.y2s;            %factor to convert years to seconds
    start = (datetime(1)-time(1))/y2s; %startime of model run (yr)
    switch option
        case 1    %linear rate of slr
            slr = dslr*time/y2s;
        case 2    %exponential from 1900
            if isempty(dslr)
                dslr = 0.011;          %exponential rate of increase
            end  
            %exponential increase post 1900 that approximates the Defra guidance
            %slr = exp(etr*T) post 1900, where T is in years 
            %For a value of dslr0 = 1mm/yr: 
            %slr is 1mm/yr in 1900; 3mm/yr in 2000; 5mm/yr in 2050; 9mm/yr in 2100 
            yr0 = 1900;                %year in which exponential increase starts (yr)
            dslr0 = 0.001;             %linear rate of slr pre 1900
            year = datetime/y2s;       %date time in years
            if year<yr0
                slr = dslr0*time/y2s;
            else
                slr0 = 0;
                if start<=yr0
                    slr0 = dslr0*(yr0-start);
                end
                slr = slr0+dslr0/dslr*(exp(dslr*(year-yr0))-1);
                slrstart = slr0+dslr0/dslr*(exp(dslr*(start-yr0))-1); %pivot year
                slr = slr-slrstart;
            end
        case 3
            %user option
        otherwise
            slr = [];
            warndlg('Undefined option in slr_rate')
    end
end
