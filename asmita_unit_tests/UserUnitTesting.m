function UserUnitTesting()
    %function to select and call selected unit test of Asmita using
    %AsmitaTest.m and runasmitamodel.m. runaasmitamodel calls Asmita in
    %silent mode and runs model for selected case and compares the results
    %with the file in muiASM TestData.

    %setting App paths - not needed
    appinfo = matlab.apputil.getInstalledAppInfo;
    idx = find(strcmp({appinfo.name},'Asmita'));
    aspath = appinfo(idx(1)).location;
    astest = [aspath,[filesep,'Asmita',filesep,'asmita_unit_tests',filesep,'AsmitaTest.m']];
    
    listtext = {'H1EMTest','HumberTest','AmelanderTest','VeniceTest','SevernTest',...
                    'InletTest','SouthamptonTest','YangtzeTest'};

    [selection, ok] = listdlg('PromptString','Select a test case:',...
                            'SelectionMode','single',...
                            'ListString',listtext);
    if ok==0, return; end %user selected cancel

    import matlab.unittest.selectors.HasParameter
    s1 = matlab.unittest.TestSuite.fromFile(astest,...
        HasParameter('Name',listtext{selection}));

    s1.run
end