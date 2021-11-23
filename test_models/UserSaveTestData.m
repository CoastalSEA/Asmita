function UserSaveTestData
%function to save a test data set for use in AsmitaTest unit testing
%This replicates the way the model is run in AsmitaTest & runasmitamodel
%Start with completed model setup and all scenarios cleared. Run is saved
%as Scenario No.1 and read in Asmita Test to check the results

%set up path to Asmita model and then call model
    aspath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\';
    userPath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\muiASM TestData\';
    testPath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\muiASM model files\';
    addpath(aspath,testPath,userPath);

    listtext = {'H1EMTest','Humber 1EM.mat';...
                'HumberTest','Humber 3EM.mat';...
                'AmelanderTest','Amelander_cfFortran.mat';...
                'VeniceTest','Venice 9EM.mat';...
                'SevernTest','Severn 16EM.mat';...
                'InletTest','Inlet 7EMts.mat';...
                'SouthamptonTest','Southampton 4EM.mat';...
                'YangtzeTest','Yangtze 30EM allInts.mat'};
                        
    [selection, ok] = listdlg('PromptString','Select a test case:',...
                            'SelectionMode','single',...
                            'ListString',listtext(:,1));
    if ok==0, return; end %user selected cancel
    
    testModel = listtext{selection,2};
    
    %create instance of asmita and load definition file
    mobj = Asmita(true);
    loadAsmitaModel(mobj,testPath,testModel);
    
    %run model and save results    
    AsmitaModel.runModel(mobj);
    
    sobj = mobj;
    sfile = sprintf('%sData',listtext{selection,1});
    save([userPath,sfile],'sobj');
    clear sobj

    %tidy up
    delete(mobj.mUI.Figure);
    delete(mobj);    %delete the class object
    clear mobj
    rmpath(userPath,testPath)
    clear aspath userPath testPath
end