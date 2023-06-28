function UserUnitTesting()
    %function to select and call selected unit test of Asmita using
    %AsmitaTest.m and runasmitamodel.m. runaasmitamodel calls Asmita in
    %silent mode and runs model for selected case and compares the results
    %with the file in muiASM TestData.

    %set up path to Asmita and then call model
    appinfo = matlab.apputil.getInstalledAppInfo;
    idx = find(strcmp({appinfo.name},'Asmita'));
    aspath = appinfo(idx(1)).location;
    if isempty(idx)                         %asmita path not found
        warndlg('Asmita App not installed'); return;
    end
    testpath = [aspath,[filesep,'Asmita',filesep,'test_models']];
    datapath = [testpath,[filesep,'muiASM TestData']];
    modelpath = [testpath,[filesep,'muiASM model files']];
    addpath(aspath,testpath,datapath,modelpath);
    
    listtext = {'H1EMTest','HumberTest','AmelanderTest','VeniceTest','SevernTest',...
                    'InletTest','SouthamptonTest','YangtzeTest'};

    [selection, ok] = listdlg('PromptString','Select a test case:',...
                            'SelectionMode','single',...
                            'ListString',listtext);
    if ok==0, return; end %user selected cancel

    import matlab.unittest.selectors.HasParameter
    s1 = matlab.unittest.TestSuite.fromFile('AsmitaTest.m',...
        HasParameter('Name',listtext{selection}));

    s1.run
    rmpath(aspath,testpath,datapath,modelpath)
    clear aspath testpath datapath modelpath
end