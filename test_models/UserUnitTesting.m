function UserUnitTesting()
    %function to select and call selected unit test of asmita

    %set up path to Asmita and then call model
    aspath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\';
    testpath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models';
    datapath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\AsmitaOO Data files';
%     datapath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\muiASM TestData';
    modelpath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\muiASM model files';
    addpath(aspath,testpath,datapath,modelpath);
    

    listtext = {'HumberTest','AmelanderTest','VeniceTest','SevernTest',...
                    'PaghamTest','InletTest'};

    [selection, ok] = listdlg('PromptString','Select a test case:',...
                            'SelectionMode','single',...
                            'ListString',listtext);
    if ok==0, return; end %user selected cancel

    import matlab.unittest.selectors.HasParameter
    s1 = matlab.unittest.TestSuite.fromFile('AsmitaTest.m',...
        HasParameter('Name',listtext{selection}));
    
%     s1 = matlab.unittest.TestSuite.fromFile('AsmitaCFmatfiles.m',...
%         HasParameter('Name',listtext{selection}));
    
    s1.run
    rmpath(aspath,testpath,datapath,modelpath)
    clear aspath testpath datapath modelpath
end
            