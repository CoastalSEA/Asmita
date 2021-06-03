function [header,data] = readinputdata()
%
%-------function help------------------------------------------------------
% NAME
%   readinputdata.m
% PURPOSE
%   read advection timeseries data from file
%   read format is file specific but can handle multiple columns
%   with space or tab separators
% USAGE
%   [header,data] = readinputdata()h
% OUTPUT
%   header - file header text 
%   data - cell array for each column of data
% SEE ALSO
%   used in Advection, Element 
%
% Author: Ian Townend
% CoastalSEA (c)Apr 2021
%--------------------------------------------------------------------------
%
    userprompt = 'Select Data file>';
    [fname, path] = uigetfile('*.txt',userprompt);            
    if fname==0   %user aborted 
        header = ''; data = [];
        return; 
    end   
    fid = fopen([path,fname], 'r');
    if fid<0
        errordlg('Could not open file for reading','File write error','modal')
        header = ''; data = [];
        return;
    end
    %find number of columns in file
    header = fgets(fid); %read first line of file       
    header = strtrim(header);
    ncol = length(strsplit(header,{'\t',' '}));
    if ncol<2
        warndlg('File should be space or tab separated columns');
        return;
    end
    dataSpec = repmat('%f ',1,ncol);
    %read numeric data                      
    data = textscan(fid,dataSpec);
    fclose(fid);
end