function UserUnitTesting()
    %function to select and call selected unit test of asmita

    %set up path to Asmita and then call model
    testpath = pwd;
    modelpath = '../';
    addpath(testpath,modelpath);
    clear testpath modelpath

    listtext = {'HumberTest','AmelanderTest','VeniceTest','SevernTest',...
                    'PaghamTest','InletTest'};

    [selection, ok] = listdlg('PromptString','Select a test case:',...
                            'SelectionMode','single',...
                            'ListString',listtext);
    if ok==0, return; end %user selected cancel

    import matlab.unittest.selectors.HasParameter
%     s1 = matlab.unittest.TestSuite.fromFile('AsmitaTest.m',...
%         HasParameter('Name',listtext{selection}));
    
    s1 = matlab.unittest.TestSuite.fromFile('AsmitaCFmatfiles.m',...
        HasParameter('Name',listtext{selection}));
    
    s1.run
end
            