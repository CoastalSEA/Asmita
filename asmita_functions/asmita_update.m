function asmita_update(obj,oldV,newV)
%
%-------header-------------------------------------------------------------
% NAME
%   asmita_update.m 
% PURPOSE
%   update saved models to newer versions of Asmita
% USAGE
%   asmita_update(oldV,newV) 
% INPUTS
%   obj - instance of model
%   oldV - old version number as a character string
%   newV - new version number as a character string
% RESULTS
%   saved model updated to new version. If this is called from Asmita this
%   will be the version that is being run at the time.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2023 
%--------------------------------------------------------------------------
%
    if strcmp(oldV,'3.30') && strcmp(newV,'3.40')
        update_v33_to_v34(obj);
    elseif strcmp(oldV,'3.30') && strcmp(newV,'3.50')
        update_v33_to_v34(obj);
        update_v34_to_v35(obj);
    elseif strcmp(oldV,'3.40') && strcmp(newV,'3.50')
        update_v34_to_v35(obj);    
    elseif strcmp(oldV,'3.30') && strcmp(newV,'4.00')
        update_v33_to_v34(obj);
        update_v34_to_v40(obj);
    elseif strcmp(oldV,'3.40') && strcmp(newV,'4.00')
        update_v34_to_v40(obj);
    elseif strcmp(oldV,'3.50') && strcmp(newV,'4.00')
        update_v34_to_v40(obj);   
    elseif strcmp(oldV,'4.00') && strcmp(newV,'3.50')    
        update_v40_to_v35(obj);
    else
        warndlg(sprintf('No update for version %s to version %s', oldV,newV))
    end
end

