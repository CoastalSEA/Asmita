function actSolution = runasmitamodel(testModel)
    %initialise Asmita, load selected test model and get results            
    mobj = Asmita(true);
    testPath = [pwd,'\muiASM model files\'];
    loadAsmitaModel(mobj,testPath,testModel);
    %run model
    AsmitaModel.runModel(mobj);
    useCase = length(mobj.Cases.DataSets.AsmitaModel);
    dst = getDataset(mobj.Cases,useCase,1);
    %first 6 properties for all time steps, all elements 
    actSolution = dst.DataTable{:,1:6};
    %need to remove totals column from each variable
    ncol = size(actSolution,2);
    ntot = ncol/6;
    nint = ntot:ntot:ncol;
    actSolution(:,nint) = [];
    
    %tidy up
    delete(mobj.mUI.Figure);
    delete(mobj);    %delete the class object 
    clear mobj
end