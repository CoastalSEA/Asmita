function [actSolution,varlist,ele] = runasmitamodel(testModel,isoldmodel)
    %initialise Asmita, load selected test model and get results            
    mobj = Asmita(true);
    testPath = [pwd,'\muiASM model files\'];
    loadAsmitaModel(mobj,testPath,testModel);
    %run model
    AsmitaModel.runModel(mobj);
    useCase = length(mobj.Cases.DataSets.AsmitaModel);
    dst = getDataset(mobj.Cases,useCase,1);    
    %need to remove totals column from each variable
    if isoldmodel
        %first 3 properties for all time steps, all elements 
        actSolution = dst.DataTable{:,1:3};
        ncol = size(actSolution,2);
        ntot = ncol/3;
        nint = ntot:ntot:ncol;
        actSolution(:,nint) = [];
    else
        actSolution = dst.DataTable{:,:};
    end
    varlist = dst.VariableNames;
    ele = dst.Dimensions.EleName;
    %tidy up
    delete(mobj.mUI.Figure);
    delete(mobj);    %delete the class object 
    clear mobj
end