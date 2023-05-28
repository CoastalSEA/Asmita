%example script to call several models and use the output

%locate Asmita App path if using App
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'Asmita'));
aspath = appinfo(idx(1)).location;
%selection from test_models. change to working folder to access other models
testpath = [aspath,[filesep,'Asmita',filesep,'asmita_unit_tests',filesep]];
modelPath = [testpath,[filesep,'muiASM model files',filesep]]; 
%or to run from Asmita folder
% modelPath = [pwd,'\asmita_unit_tests\muiASM model files\']; 

workfolder = pwd;
cd(modelPath)  %model to folder with files

userprompt = 'Select model file(s)>';
[fname, ~]=uigetfile('*.mat',userprompt,'MultiSelect','on');
if isequal(fname,0), return; end

cd(workfolder) %restore working folder

nrec = length(fname);                    %number of models
%loop over each model file, run model and create array of model instances
for i=1:nrec
    modelName = fname{i};
    mobj(i) = asmita_command_line_run(modelPath,modelName); %#ok<SAGROW> can't intialise because Matlab calls Asmita without true argument to fill array
end

%now have an array of model instances to do something with.----------------
hf = figure; ax = axes(hf);
legtext = cell(1,nrec);
hold on
for j=1:nrec
    %get results of run in a dstable
    % mobj(j).Cases is the jth struct of the model output cases
    useCase = length(mobj(j).Cases.DataSets.AsmitaModel);
    dst = getDataset(mobj(j).Cases,useCase,1);  %use case just addded
    vm = dst.MovingVolume(:,end);               %moving volumes
    time =dst.RowNames;
    plot(ax,time,vm)
    legtext{j} = ['Moving ',fname{j}];
end
hold off
xlabel('Time')
ylabel('Volume')
legend(legtext{:},'Location','southwest') 

%tidy up-------------------------------------------------------------------
clear modelName modelPath workfolder userprompt       %setup variables
clear anobj fname i
clear hf ax dst useCase time vm j legtext nrec        %plot variables
delete(mobj);    %delete the class object
clear mobj