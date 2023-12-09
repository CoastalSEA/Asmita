function setMultiEle(src,~) 
%
%-------function help------------------------------------------------------
% NAME
%   setMultiEle.m
% PURPOSE
%   callback function for button to set multi-element plot selection
% USAGE
%   callback function: @(src,evt)setMultiEle(src,evt)
% INPUT
%   src - handle to calling object (eg graphical button)
% OUTPUT
%   change the String, UserData and Tooltip of src object
% SEE ALSO
%   used in ASM_PlotsUI derived from muiDataUI abstract class
%
% Author: Ian Townend
% CoastalSEA (c)June 2021
%--------------------------------------------------------------------------
%
    %
    if strcmp(src.String,'AE')
        src.String = 'SE';
        src.UserData = 0;
        src. Tooltip = 'Switch to plot all elements';
    elseif strcmp(src.String,'SE')
        src.String = 'AE';
        src.UserData = 1;
        src. Tooltip = 'Switch to plot single element';
    end
end
