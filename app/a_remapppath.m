%script to remove path for use in app
apppath = pwd;
appdoc = [apppath,'/doc'];
apphelp = [apppath,'/help'];
appeg = [apppath,'/example'];
appfncs = [apppath,'/asmita_functions'];
applib1 = '../../muiAppLib/muiAppEstuaryFcns'; 

rmpath(apppath,appdoc,genpath(apphelp), appeg, genpath(appfncs), applib1)

clear  apppath appdoc apphelp appeg appfncs applib1
