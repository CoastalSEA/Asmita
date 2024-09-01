%example script to run Asmita from command prompt

%locate Asmita App path if using App
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'Asmita'));
aspath = appinfo(idx(1)).location;
testpath = [aspath,[filesep,'Asmita',filesep,'asmita_unit_tests',filesep]];
modelPath = [testpath,[filesep,'muiASM model files',filesep]];
%or to run from Asmita folder
% modelPath = [pwd,'\asmita_unit_tests\muiASM model files\']; 

%create instance of Asmita and load existing model
mobj = Asmita(true);                       %true runs Asmita in silent mode
modelName = 'Humber 3EM.mat';
loadAsmitaModel(mobj,modelPath,modelName);

%Change a property---------------------------------------------------------
% requires some knowledge of the model classes and properties - see documentation).
% mobj.Inputs is a struct of the classes used to setup the model

% e.g. change global equilbirum concentration
mobj.Inputs.Estuary.EqRhoCoarse = 0.35;
mobj.Inputs.Estuary.EqRhoFine = 0;

% e.g. change transport coefficients and vertical exchange for all intertidal elements
eleobj = mobj.Inputs.Element;                    %element object
eletype = getEleProp(eleobj,'EleType');          %element types
idx = matches(eletype,'Tidalflat');              %index to all tidal flat elements
eleobj(idx).VerticalExchange = 0.01;             %change vertical exchange
eleobj(idx).TransportCoeff = 5;                  %change transport coefficient
setClassObj(mobj,'Inputs','Element',eleobj);     %reassign to model object

% e.g. change run time step
mobj.Inputs.RunProperties.TimeStep = 0.25;       %double length of time step
mobj.Inputs.RunProperties.NumSteps = 1280;       %adjust number of steps
mobj.Inputs.RunProperties.OutInterval = 4;       %adjust output interval

% e.g. change river input sediment discharge and sediment load
% case is for a singlevalue, not a timeseries, and only a single river
% input defined. If multipe rivers need rivobj is an array of River
% instances and need to select and amend each instance. New timeseries
% require RiverTSC to be modified (tsRiverFlow and tsRiverConc are transient)
mobj.Inputs.River.RiverFlow = 500;
mobj.Inputs.River.RiverRho = 0.02;

%run model-----------------------------------------------------------------
AsmitaModel.runModel(mobj);                %all paths must be set (should be done if using App)

%get results of run in a dstable-------------------------------------------
% mobj.Cases is a struct of the model output cases
useCase = length(mobj.Cases.DataSets.AsmitaModel);
dst = getDataset(mobj.Cases,useCase,1);  

%do something with the results---------------------------------------------
% e.g. plot 'All' elements volumes for Humber 3EM
vm = dst.MovingVolume(:,4);      
vf = dst.FixedVolume(:,4);
ve = dst.EqVolume(:,4);
time =dst.RowNames;

hf = figure; ax = axes(hf);
plot(ax,time,vm)
hold on
plot(ax,time,vf)
plot(ax,time,ve)
hold off
xlabel('Time')
ylabel('Volume')
legend('Moving','Fixed','Equilibirum','Location','southwest')

%tidy up-------------------------------------------------------------------
clear modelName modelPath                        %setup variables
clear hf ax dst useCase time ve vm vf            %plot variables
clear eleobj eletype idx                         %change variables
delete(mobj);    %delete the class object
clear mobj
