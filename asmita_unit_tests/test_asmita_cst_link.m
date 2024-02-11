function test_asmita_cst_link(isApp)
    %function to check the running of the CSTmodel from within Asmita
    if nargin<1
        isApp = false;
    end

    if isApp
        %locate Asmita App path if using App
        appinfo = matlab.apputil.getInstalledAppInfo;
        idx = find(strcmp({appinfo.name},'Asmita'));
        if isempty(idx), warndlg('App not found'); return; end
        aspath = appinfo(idx(1)).location;
        testpath = [aspath,[filesep,'Asmita',filesep,'asmita_unit_tests',filesep]];
        modelPath = [testpath,[filesep,'muiASM model files',filesep]];
    else
        %or to run from Asmita root folder
        aspath = pwd;
        idx = regexp(pwd,'Asmita');
        modelPath = [aspath(1:idx+5),'\asmita_unit_tests\muiASM model files\']; 
    end
    
    %create instance of Asmita and load existing model
    obj = Asmita(true);                       %true runs Asmita in silent mode
    modelName = 'Yangtze 30EM allInts.mat';
    loadAsmitaModel(obj,modelPath,modelName);

    msgtxt = 'Hydraulic properties have not been defined';
    cobj = getClassObj(obj,'Inputs','CSThydraulics',msgtxt);
    if isempty(cobj), return; end
    runModel(cobj,obj);
    tabPlot(cobj);
end
