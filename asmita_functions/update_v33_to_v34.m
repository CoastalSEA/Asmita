function update_v33_to_v34(obj)
%
%-------header-------------------------------------------------------------
% NAME
%   update_v34_to_v34.m 
% PURPOSE
%   update saved models from v3.3 to v3.4.
% USAGE
%   update_v33_to_v34(obj)
% INPUTS
%   obj - instance of model
% RESULTS
%   saved model updated from v3.3 to v3.4.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
% To use from command line, open ASMITA using:
% >>sm=Asmita;     and then call
% >>update_v33_to_v34(sm)
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2023 
%--------------------------------------------------------------------------
%

    %modifies Saltmarsh properties and PropertyLabels
    if isfield(obj.Inputs,'Saltmarsh') && ~isprop(obj.Inputs.Saltmarsh,'EdgeErosion')
        obj.Inputs.Saltmarsh.PropertyLabels = {'Number of saltmarsh species',...
                                  'Minimum depth (m)','Maximum depth (m)',...
                                  'Maximum biomass (kg/m2)',...
                                  'Species productivity (m2/kg/yr)',...
                                  'Settling coefficient, alpha (1/s)',...
                                  'Settling exponent, beta (-)',...
                                  'Edge erosion rate (m/yr)',...
                                  'Include marsh flat erosion (0/1)'};
        
        obj.Inputs.Saltmarsh.EdgeErosion = 0;
        obj.Inputs.Saltmarsh.FlatErosion = false;
    end
    
    %modifies any defined Interventions to include non-erodibility (true is so)
    if isfield(obj.Inputs,'Interventions') && ~isprop(obj.Inputs.Interventions,'isNonErodible')
        lobj = getClassObj(obj,'Inputs','Interventions');
        nint = length(lobj);   %number of element interventions defined
        if nint>0 && iscell(lobj(1).Year)
            for i=1:nint
                nyears = length(lobj(i).Year{1});
                lobj(i).isNonErodible = {zeros(nyears,1)}; %as a cell array
            end    
            setClassObj(obj,'Inputs','Interventions',lobj);
        end    
        clear lobj nint nyears i
    end
    
    %change Interventions input from cell array to numerical array
    if isfield(obj.Inputs,'Interventions')
        intobj = getClassObj(obj,'Inputs','Interventions');
        nint = length(intobj);
        if nint>0 && iscell(intobj(1).Year)
            for i=1:nint
                obj.Inputs.Interventions(i).Year = intobj(i).Year{1};
                obj.Inputs.Interventions(i).VolumeChange = intobj(i).VolumeChange{1};
                obj.Inputs.Interventions(i).SurfaceAreaChange = intobj(i).SurfaceAreaChange{1};
                obj.Inputs.Interventions(i).isNonErodible = intobj(i).isNonErodible{1};
            end
        end
        clear intobj nint i
    end
    
    %update RunParams if Cases have been run
    if ~isempty(obj.Cases) && isfield(obj.Cases.DataSets,'AsmitaModel')
        ncase = length(obj.Cases.DataSets.AsmitaModel);
        for i=1:ncase
            if isfield(obj.Cases.DataSets.AsmitaModel(i).RunParam,'Interventions')
                localObj = obj.Cases.DataSets.AsmitaModel(i).RunParam.Interventions;
                if ~isempty(localObj) && iscell(localObj(1).Year)
                    nint = length(obj.Cases.DataSets.AsmitaModel(i).RunParam.Interventions);
                    for j=1:nint
                        %add isNonErodible with correct dimensions
                        nyears = length(localObj(j).Year{1});
                        localObj(j).isNonErodible = {zeros(nyears,1)}; %as a cell array
                        %change from cell to numerical array
                        localObj(j).Year = localObj(j).Year{1};
                        localObj(j).VolumeChange = localObj(j).VolumeChange{1};
                        localObj(j).SurfaceAreaChange = localObj(j).SurfaceAreaChange{1};
                        localObj(j).isNonErodible = localObj(j).isNonErodible{1};
                    end   
                end
                obj.Cases.DataSets.AsmitaModel(i).RunParam.Interventions = localObj;
            end
        end
        clear ncase localObj nint nyears i j
    end
    
    getdialog('Project updated from v3.3 to v3.4')
end