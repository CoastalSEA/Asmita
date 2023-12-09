classdef RunProperties < muiPropertyUI              
%
%-------class help---------------------------------------------------------
% NAME
%   RunProperties.m
% PURPOSE
%   Class for to set run time properties - used in ASMITA and ChannelForm
% USAGE
%   obj = RunProperties.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {' Time Step (years)',...
                          ' Number of Time Steps',...
                          ' Output Interval (No. of time steps)', ...
                          ' Start Year',...
                          ' Run time plot'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        TimeStep = 1         %timestep in years (can be decimal part)
        NumSteps = 1         %number of time steps to be used in run
        OutInterval = 1      %number of step intervalss for output sub-sampling
        StartYear = 1900     %start year for model 
        isRunPlot = false    %option to have a plot of volume prism and area at run time
    end    

%%   
    methods (Access=protected)
        function obj = RunProperties(mobj)         
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'RunProperties';           
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = RunProperties(mobj);             
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
    end
end