function UserSaveTestData(aspath)
    %function to save a test data set for use in AsmitaTest unit testing
    %This replicates the way the model is run in AsmitaTest & runasmitamodel
    %Start with completed model setup and all scenarios cleared. Run is saved
    %as Case No.1 and read in Asmita Test to check the results (NB. Can
    %also do this using the UI).
   
    %setting App paths
    if nargin<1
        appinfo = matlab.apputil.getInstalledAppInfo;
        idx = find(strcmp({appinfo.name},'Asmita'));
        aspath = appinfo(idx(1)).location;
        if isempty(idx)                         %asmita path not found
            warndlg('Asmita App not installed'); return;
        end
    end
    testpath = [aspath,[filesep,'Asmita',filesep,'asmita_unit_tests']];
    userPath = [testpath,[filesep,'muiASM TestData',filesep]];
    testPath = [testpath,[filesep,'muiASM model files',filesep]];

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
    delete(mobj);    %delete the class object
    clear mobj

    getdialog(sprintf('File: %s has been saved',sfile))
end