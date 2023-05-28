function asmitaOO_structfile()
%
%-------function help------------------------------------------------------
% NAME
%   asmitaOO_structfile.m
% PURPOSE
%   save AsmitaOO projects as a struct called 'inp' htat uses the same
%   class handles and property names as sobj
% USAGE
%   asmitaOO_structfile()
% NOTES
%   user is prompted to select an input mat file and the converted project
%   is saved to a new mat file containing the struct called 'inp'
%   in order to read the correct properties for each class the paths for
%   source code of AsmitaOO ensure that only the AsmitaOO path is
%   accessible (ie the path to the mui version of Asmita is not set)
% SEE ALSO
%   asm_oo2mui.m AsmitaOO and muiApp/Asmita
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    %load the AsmitaOO version of the model    
    asmOOpath = 'D:\Work\Tools\MATLAB\AsmitaOO';
    addpath(asmOOpath)

    [sfile,spath] = uigetfile({'*.mat','MAT-files (*.mat)'},'Open AsmitaOO mat file');
    if sfile==0, return; end
    
    load([spath,sfile],'sobj');
    fnames = fieldnames(sobj);  
    
    for ii=1:length(fnames)
        lobj = sobj.(fnames{ii});
        if ~isempty(lobj)
            propnames = getPropertyNames(lobj);
            for jj=1:length(lobj)
                for kk=1:length(propnames)                
                    inp.(fnames{ii})(jj).(propnames{kk}) = lobj(jj).(propnames{kk});
                end
            end
        end
    end
    
    rmpath(asmOOpath)
    clear sobj
    
    sfile = sprintf('inp%s',sfile);
    save(['D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\AsmitaOO inp files\',sfile],'inp');
    msgbox(sprintf('File saved as %s',sfile));
end
%%
function propnames = getPropertyNames(obj)
    %for class obj get the property names that are not Transient 
    %or Hidden
    mc = metaclass(obj);
    mp = mc.PropertyList;
%     ms = mp(1).DefiningClass.SuperclassList;
%     scnames = getSuperclassNames(obj,ms);
    count=1;
    propnames = {};
    for k=1:length(mp)
        %remove hidden, transient and constant properties
        idx = mp(k).Transient + mp(k).Hidden + mp(k).Constant +...
                                               + mp(k).Dependent;
                                           
%         idx = mp(k).Constant + mp(k).Dependent;                              
%         %remove superclass properties
%         idx = idx + any(strcmp(scnames, mp(k).Name)); 
        if idx<1
            propnames{count,1} = mp(k).Name;
            count = count+1;
        end
    end
    %
end