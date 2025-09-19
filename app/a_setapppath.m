%script to set path for use in app
%run from the ../Asmita/app folder with supporting muiAppLib, dstoolbox and
%muitoolbox branching from the same parent folder (eg muiModels)
apppath = pwd;
appdoc = [apppath,'/doc'];
apphelp = [apppath,'/help/html'];
appeg = [apppath,'/example'];
appfncs = [apppath,'/asmita_functions'];
applib1 = '../../muiAppLib/muiAppEstuaryFcns'; 

addpath(apppath, appdoc, genpath(apphelp), appeg, genpath(appfncs), applib1)

clear apppath appdoc apphelp appeg appfncs applib1 

%if toolboxes not installed
isdstbx =  which('dstable.m');
if isempty(isdstbx)
    dstbx = 'D:\Matlab Code\muiModels\dstoolbox\toolbox';
    addpath(genpath(dstbx));
    clear dstbx 
end

ismuitbx =  which('muiModelUI.m');
if isempty(ismuitbx)
    muitbx = 'D:\Matlab Code\muiModels\muitoolbox\toolbox';
    addpath(genpath(muitbx))
    clear muitbx
end
clear isdstbx ismuitbx