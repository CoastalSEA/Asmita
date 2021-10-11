classdef RunConditions < muiPropertyUI                
%
%-------class help---------------------------------------------------------
% NAME
%   RunConditions.m
% PURPOSE
%   Class for defining the ASMITA run time conditions
% USAGE
%   obj = RunConditions.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {' Scale to initial values ', ...
                          ' Include interventions ', ...
                          ' Include river advection ', ...
                          ' Include river flow offset ', ...
                          ' Include drift advection ', ...
                          ' Include drift rate offset ', ...
                          ' Include saltmarsh biomass ',...
                          ' Include tidal pumping ',...
                          ' Include dynamic hydraulics'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        DriftTSC              %timeseries collection of drift rate (m3/year)
        IncDynamicElements = false %not active but needed because used in code        
    end
    
    properties
        ScaleValues = false     %flag to include scaling of initial values
        IncInterventions = false %flag to include user defined interventions
        IncRiver = false        %flag to include river advections
        RiverOffset = false     %flag to correct equilbirium values for river flow
        IncDrift = false        %flag to include drift advections
        DriftOffset = false     %flag to correct equilbirium values for drift
        IncSaltmarsh = false    %flag to include saltmarsh biological productivity
        IncTidalPumping = false %flag to include tidal pumping
        IncDynHydraulics = false%flag to include hydraulics data from CSTmodel
%         IncDynamicElements = false %flag to enable elements to change type and connectivity        
    end   
    
    properties (Transient)
        Adv2Inc                 %advections to include in model run
        Adv2Offset              %advections to offset in initial equilibrium              
    end
%%   
    methods (Access=protected)
        function obj = RunConditions(mobj)   
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
            classname = 'RunConditions';  
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = RunConditions(mobj);             
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end  
%%
        function setAdvectionOffset(mobj)
            %find whether an offset is to be included
            obj = getClassObj(mobj,'Inputs','RunConditions');
            if obj.IncRiver && obj.IncDrift %both included
                obj.Adv2Inc = 'flow+drift';           %include both;
                if obj.RiverOffset && obj.DriftOffset
                    obj.Adv2Offset = 'flow+drift';        %include both offsets
                elseif obj.RiverOffset
                    obj.Adv2Offset = 'flow';              %include river offset
                elseif obj.DriftOffset
                    obj.Adv2Offset = 'drift';             %include drift offset
                else
                    obj.Adv2Offset = 'none';
                end         
            elseif obj.IncRiver                %only river included
                obj.Adv2Inc = 'flow';                 %include river
                if obj.RiverOffset
                    obj.Adv2Offset = 'flow';              %include river offset
                else
                    obj.Adv2Offset = 'none';
                end
            elseif obj.IncDrift                %only drift included
                obj.Adv2Inc = 'drift';                %include drift
                if obj.DriftOffset
                    obj.Adv2Offset = 'drift';             %include drift offset
                else
                    obj.Adv2Offset = 'none';    
                end
            else
                obj.Adv2Inc = 'none';
                obj.Adv2Offset = 'none';
            end
            setClassObj(mobj,'Inputs','RunConditions',obj);
        end
    end
%%        
        %add other functions to operate on properties as required   
end