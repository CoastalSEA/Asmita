function update_v34_to_v35(obj)
%
%-------header-------------------------------------------------------------
% NAME
%   update_v34_to_v35.m 
% PURPOSE
%   update saved models from v3.4 to v3.5.
% USAGE
%   update_v34_to_v35(obj)
% INPUTS
%   obj - instance of model
% RESULTS
%   saved model updated from v3.4 to v3.5.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
% To use from command line, open ASMITA using:
% >>sm=Asmita;     and then call
% >>update_v34_to_v35(sm)
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2023 
%--------------------------------------------------------------------------
%

    %add run time plot option to RunProperties and sadjust tab position
    if isfield(obj.Inputs,'RunProperties') && ...
                            ~isfield(obj.Inputs.RunProperties,'isRunPlot')
        obj.Inputs.RunProperties.PropertyLabels = {' Time Step (years)',...
                          ' Number of Time Steps',...
                          ' Output Interval (No. of time steps)', ...
                          ' Start Year',...
                          ' Run time plot (0/1)'};

        obj.Inputs.RunProperties.isRunPlot = false;
        obj.Inputs.RunProperties.TabDisplay.Position(1) = 0.94;        
    end

    %adjust tab position for EqCoeffParams
    if isfield(obj.Inputs,'EqCoeffParams')
        obj.Inputs.EqCoeffParams.TabDisplay.Position(1) = 0.94;
    end
    
     getdialog('Project updated from v3.4 to v3.5')
end