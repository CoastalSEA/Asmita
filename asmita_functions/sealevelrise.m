function [slr,dslr] = sealevelrise(yrs,pivotyear,dslrvec,option)
%
%-------function help------------------------------------------------------
% NAME
%   sealevelrise.m
% PURPOSE
%   Function to compute sea level rise 
% USAGE
%   slr = sealevelrise(time,datetime,dslr,option)
% INPUTS
%   yrs - vector of time in years
%   pivotyear - pivot year for zero sea level rise (ie negative for years 
%               before pivot year) - note this is not the same as yr0
%   dslrvec - vector containing rate parameters
%              linear - rate of sea level rise (m/yr)
%              exponential - [dslr, yr0, dslr0], where:
%                  dsr - exponential rate increase (m/yr), 
%                  yr0 - start of exponential (year), 
%                  dslr0 - linear rate pre start of exponential
%              if dslrvec is empty default values used [0.011,1900,0.001]             
%   option - determine how slr is calculated
%              1 - linear slr
%              2 - exponential from 1900
%              3, etc - user can add addtional options as required
% OUTPUTS
%   slr = magnitude of sea level rise from start time to time now (m)
%   dslr - rate of sea level rise at time t (m/yr)
%          only changes if rate is not linear, otherwise = dslrvec
% EXAMPLE
%   yrs = 1800:1:2100;
%   [slr,dslr] = sealevelrise(yrs,2000,[],2); %exponential slr from 1900,
%                                             %with zero in the year 2000 
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    switch option
        case 1    %linear rate of slr
            elapsedyears = yrs-pivotyear;
            slr = dslrvec*elapsedyears;
            dslr = dslrvec;
        case 2    %exponential rate of slr
            if isempty(dslrvec) || length(dslrvec)<3
                %exponential increase post 1900 that approximates the Defra guidance
                %slr = exp(etr*T) post 1900, where T is in years 
                %For a value of dslr0 = 1mm/yr: 
                %  slr is 1mm/yr in 1900; 3mm/yr in 2000; 5mm/yr in 2050; 
                %  9mm/yr in 2100 
                expslr = 0.011;      %exponential rate of increase
                yr0 = 1900;          %year in which exponential increase starts (yr)
                dslr0 = 0.001;       %linear rate of slr pre 1900
            else
                expslr = dslrvec(1); %exponential rate of increase
                yr0 = dslrvec(2);    %year in which exponential increase starts (yr)
                dslr0 = dslrvec(3);  %linear rate of slr pre 1900
            end
            
            %compute any linear component of slr
            slr = zeros(size(yrs));            
            idx = yrs>yr0;   %index of records before switch to exponential
            slr(~idx) = dslr0*(yrs(~idx)-yrs(1)); %linear slr

            
            %compute the exponential
            slr0 = dslr0*(yr0-yrs(1));            %slr at yr0
            slr(idx) = slr0+dslr0/expslr*(exp(expslr*(yrs(idx)-yr0))-1);
            
            %compute the rate of slr
            dslr = [dslr0,diff(slr)./diff(yrs)];
            
            %offset the curve so slr=0 in the pivotyear
            if pivotyear<yr0
                slroffset = dslr0*(pivotyear-yrs(1)); 
            else
                slroffset = slr0+dslr0/expslr*(exp(expslr*(pivotyear-yr0))-1);   
            end
            slr = slr-slroffset;            
        case 3
            %user option
        otherwise
            slr = []; dslr = [];
            warndlg('Undefined option in sealevelrise')            
    end
end
            