function mobj = asmita_command_line_run(modelPath,modelName)
%
%-------function help------------------------------------------------------
% NAME
%   asmita_command_line_run.m
% PURPOSE
%   run the Asmita model without opening the UI and without any run time
%   messages
% USAGE
%   mobj = asmita_command_line_run(modelPath,modelName)
% INPUTS
%   modelPath - path to folder that contains the model to be run
%   modelName - name of model mat file to use
% OUTPUTS
%   mobj - model instance with no UI
% SEE ALSO
%   asmita command line script.m
%
% Author: Ian Townend
% CoastalSEA (c) June 2023
%--------------------------------------------------------------------------
%

    %create instance of Asmita and load existing model
    mobj = Asmita(true);         %true runs Asmita in silent mode
    loadAsmitaModel(mobj,modelPath,modelName);   
    %run model-------------------------------------------------------------
    AsmitaModel.runModel(mobj);  %all paths must be set (should be done if using App)
end