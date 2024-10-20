function asmitaOO_structfile()
%
%-------function help------------------------------------------------------
% NAME
%   asmitaOO_structfile.m
% PURPOSE
%   oo2inp - save AsmitaOO projects as a struct called 'inp' that uses the 
%   same class handles and property names as sobj
% USAGE
%   asmitaOO_structfile()
% NOTES
%   user is prompted to select an asmitaOO mat file and the converted project
%   is saved to a new 'inp' mat file, containing the struct called 'inp'
%   in order to read the correct properties for each class the paths for
%   source code of AsmitaOO ensure that only the AsmitaOO path is
%   accessible (ie the path to the mui version of Asmita is not set)
% SEE ALSO
%   asm_oo2mui.m, AsmitaOO and muiApp/Asmita
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021 and May 2023
%--------------------------------------------------------------------------
%
    %check whether App path is set and if so warn and quit 
    appinfo = matlab.apputil.getInstalledAppInfo;
    idx = find(strcmp({appinfo.name},'Asmita'));
    s = pathsep;
    pathStr = [s, path, s];
    if ~isempty(idx)                         %asmita path found
        aspath = appinfo(idx(1)).location;        
        onPath  = contains(pathStr, [s, aspath, s], 'IgnoreCase', ispc); %ignores case on Windows pc 
        if onPath                            %asmita path in Matlab search path
            warndlg('Remove Asmita App paths before running OO2inp script')
            return;
        end
    end

    asPath = contains(pathStr, [s,'AsmitaOO', s], 'IgnoreCase', ispc); %ignores case on Windows pc
    if ~asPath                               %set path to AsmitaOO source code
        asmOOpath = 'D:\OneDrive\Software_Tools\Tools\MATLAB\AsmitaOO';
        answer = inputdlg('AsmitaOO location:','asmOOpath',1,{asmOOpath});
        if isempty(answer), return; end
        asmOOpath = answer{1};
        addpath(asmOOpath)
    end
    [sfile,spath] = uigetfile({'*.mat','MAT-files (*.mat)'},'Open AsmitaOO mat file');
    if sfile==0, return; end
    
    %create inp struct with properties from AsmitaOO classes
    load([spath,sfile],'sobj');
    fnames = fieldnames(sobj);  
    
    inp.ProjectInfo = sobj.ProjectInfo;
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
    
    %save the resulting inp struct
    sfile = sprintf('%s_inp.mat',sfile(1:end-4));
    save([pwd,filesep,sfile],'inp');
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