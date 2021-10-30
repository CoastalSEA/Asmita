function asm = asm_oo2mui()
%
%-------function help------------------------------------------------------
% NAME
%   asm_oo2mui.m
% PURPOSE
%   convert AsmitaOO project files to the format used by the version of
%   Asmita based on the muitoolbox
% USAGE
%   asm_oo2mui()
% NOTES
%   user is prompted to select an input mat file and the converted project
%   is saved to a new mat file
%   in order to read the correct properties for each class the paths for
%   source code of AsmitaOO and the mui version of Asmita need to be edited
%   at the 'cd' commands below 
% SEE ALSO
%   AsmitaOO and muiApp/Asmita
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    [sfile,spath] = uigetfile({'*.mat','MAT-files (*.mat)'},'Open AsmitaOO mat file');
    if sfile==0, return; end

    %load the AsmitaOO version of the model    
%     asmOOpath = 'D:\Work\Tools\MATLAB\AsmitaOO';
%     cd 'D:\Work\Tools\MATLAB\AsmitaOO'
%     addpath(asmOOpath)
%     rmpath('D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita')
    
    load([spath,sfile],'inp');
%     fnames = fieldnames(sobj);  
%     
%     for ii=1:length(fnames)
%         lobj = sobj.(fnames{ii});
%         if ~isempty(lobj)
%             propnames = getPropertyNames(lobj);
%             for jj=1:length(lobj)
%                 for kk=1:length(propnames)                
%                     inp.(fnames{ii})(jj).(propnames{kk}) = lobj(jj).(propnames{kk});
%                 end
%             end
%         end
%     end
%     
%     rmpath(asmOOpath)
%     clear sobj asmOOpath
%     
%     cd 'D:\Work\Tools\MATLAB\MUImodels2'
%     muisetpaths
%     cd 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita'
    addpath('D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita')
    %open an instance of the mui version of Asmita in silent mode    
    asm = Asmita(true);  %true flag runs Asmita in silent mode 
    
    %now unpack old model and load into new model
    if ~isempty(inp.ProjectInfo)
        asm.Info.ProjectName  = inp.ProjectInfo.ProjectName;
        asm.Info.ProjectDate = inp.ProjectInfo.ProjectDate;
    end
    %assume constants are unchanged
    
    %handle each of the old classes in turn
    fnames = fieldnames(inp);    
    for kf=3:length(fnames)
        if ~isempty(inp.(fnames{kf}))
            switch fnames{kf}
                case 'h_wlv'
                    old = inp.h_wlv;
                    wlv = WaterLevels.setInput(asm,false);
                    sameprops = {'TidalPeriod','LWtoHWratio','MSL0',...
                        'SLRrate','M2amplitude','S2amplitude',...
                        'M4amplitude','M4phase'};
                    for j= 1:length(sameprops)
                        wlv.(sameprops{j}) = old.(sameprops{j});
                    end
                    wlv.TidalAmp = old.TidalRange/2;      %tidal amplitude (m)
                    wlv.NumCycles = old.RangeCycles.Num;  %number of cycles to include
                    wlv.CycleAmp = old.RangeCycles.Amp';   %vector of cycle amplitudes (m)
                    wlv.CyclePeriod = old.RangeCycles.Per';%vector of cycle periods (yrs)
                    wlv.CyclePhase = old.RangeCycles.Pha'; %vector of cycle phases (yrs)
                    asm.Inputs.WaterLevels = wlv;
                    
                case 'h_rnp'
                    %deal with RunProperties first
                    old = inp.h_rnp;
                    rnp = RunProperties.setInput(asm,false);
                    sameprops = {'TimeStep','NumSteps','OutInterval','StartYear'};
                    for j= 1:length(sameprops)
                        rnp.(sameprops{j}) = old.(sameprops{j});
                    end
                    asm.Inputs.RunProperties = rnp;
                    
                    %now update RunConditions
                    rnc = RunConditions.setInput(asm,false);
                    sameprops = {'ScaleValues','RiverOffset','DriftOffset',...
                        'IncSaltmarsh','IncTidalPumping'};
                    for j= 1:length(sameprops)
                        rnc.(sameprops{j}) = old.(sameprops{j});
                    end
                    rnc.IncInterventions = old.IncIntervention;
                    %amend to match new condition settings
                    if rnc.RiverOffset, rnc.IncRiver = true; end
                    if rnc.IncTidalPumping, rnc.IncRiver = true; end
                    if rnc.DriftOffset, rnc.IncDrift = true; end
                    asm.Inputs.RunConditions = rnc;
                    
                    %now update EqCoeffProps
                    eqp = EqCoeffProps.setInput(asm,false);
                    sameprops = {'alpha','beta','eqtype'};
                    for j= 1:length(sameprops)
                        eqp.(sameprops{j}) = old.(sameprops{j});
                    end
                    asm.Inputs.EqCoeffProps = eqp;
                    
                case 'h_est'
                    old = inp.h_est;
                    est = Estuary.setInput(asm,false);
                    sameprops = {'Dispersion','ExternalDisp','WidthELength',...
                        'WindSpeed','WindHeight'};
                    for j= 1:length(sameprops)
                        est.(sameprops{j}) = old.(sameprops{j});
                    end
                    est.EqRhoCoarse = old.EqConcCoarse*2650;
                    est.EqRhoFine = old.EqConcFine*2650;
                    asm.Inputs.Estuary = est;
                    
                case 'h_ele'
                    old = inp.h_ele;
                    sameprops = {'EleID','ReachID','EleType','EleName',...
                        'Erodible','InitialVolume','InitialSurfaceArea',...
                        'Length','VerticalExchange','TransportCoeff',...
                        'SedMobility','TidalDamping'};
                    nele = length(old);
                    Element.setElement(asm);
                    anele = getClassObj(asm,'Inputs','Element');
                    for iele=1:nele
                        ele(iele) = copy(anele); %#ok<*AGROW>
                        for j=1:length(sameprops)
                            if iscell(old(iele).(sameprops{j}))
                                ele(iele).(sameprops{j}) = old(iele).(sameprops{j}){1};
                            else
                                ele(iele).(sameprops{j}) = old(iele).(sameprops{j});
                            end
                        end
                        ele(iele).BedDensity = old(iele).BedConcentration*2650;
                    end
                    asm.Inputs.Element = ele;
                    
                case 'h_riv'
                    old = inp.h_riv;
                    sameprops = {'ChannelID','RiverFlow','RiverTSC'}; 
                    nriv = length(old);
                    River.setNewRiver(asm);
                    ariv = getClassObj(asm,'Inputs','River');
                    for iriv=1:nriv
                        riv(iriv) = copy(ariv); 
                        for j= 1:length(sameprops)                        
                            riv(iriv).(sameprops{j}) = old(iriv).(sameprops{j});
                        end
                        riv(iriv).RiverRho = old(iriv).RiverConc*2650;
                    end
                    asm.Inputs.River = riv;
                    
                case 'h_int'
                    old = inp.h_int;
                    sameprops = {'ElementID','Year','VolumeChange',...
                                 'SurfaceAreaChange'};
                    nint = length(old);
                    Interventions.setNewIntervention(asm);
                    aint = getClassObj(asm,'Inputs','Interventions');
                    for iint=1:nint
                        int(iint) = copy(aint); 
                        for j= 1:length(sameprops)
                            int(iint).(sameprops{j}) = old(iint).(sameprops{j});
                        end
                    end
                    asm.Inputs.Interventions = int;
                    
                case 'h_dft'
                    old = inp.h_dft;
                    sameprops = {'DriftEleID','DriftRate','DriftTSC'};
                    ndft = length(old);
                    Drift.setNewDrift(asm);
                    adft = getClassObj(asm,'Inputs','Drift');
                    for idft=1:ndft
                        dft(idft) = copy(adft); 
                        for j= 1:length(sameprops)
                            dft(idft).(sameprops{j}) = old(idft).(sameprops{j});
                        end
                    end
                    asm.Inputs.Drift = dft;
                    
                case 'h_sms'
                    old = inp.h_sms;
                    Saltmarsh.setInput(asm,false);
                    sms = getClassObj(asm,'Inputs','Saltmarsh');
                    sameprops = {'NumSpecies','MinSpDepth','MaxSpDepth',...
                                 'MaxBiomass','SpeciesProduct','SettlingAlpha',...
                                 'SettlingBeta','MinEdgeErosion','MaxEdgeErosion'};                        
                    for j= 1:length(sameprops)
                        sms.(sameprops{j}) = old.(sameprops{j});
                    end
                    asm.Inputs.Saltmarsh = sms;
                    
                case 'h_adv'
                    old = inp.h_adv;
                    Advection.setNewAdvection(asm);
                    adv = getClassObj(asm,'Inputs','Advection');
                    sameprops = {'RiverFlows','RiverIn','RiverOut',...
                                 'DriftFlows','DriftIn','DriftOut'};
                    for j= 1:length(sameprops)
                        adv.(sameprops{j}) = old.(sameprops{j});
                    end
                    asm.Inputs.Advection = adv;
            end
        end
    end
    
    %save new model to mat file
    asm.Info.PathName = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\test_models\muiASM model files';
    asm.Info.FileName = sfile(4:end);
    saveAfile(asm);
    closeMainFig(asm) 
    rmpath('D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita')
end


%%        
function obj = saveAfile(obj)
    %prompt user to Save As in selected folder
    pname = obj.Info.PathName;
    ispath = false;
    if ~isempty(pname)   %check if already using a path
        try
            cpath = cd(pname);
            ispath = true;
        catch
            getdialog('Path not found')
        end
    end

    [sfile,spath] = uiputfile({'*.mat','MAT-files (*.mat)'}, ...
                   'Save Asmita work file');
    if isequal(sfile,0)|| isequal(spath,0)
        return     %user selected cancel
    else
        obj.Info.PathName = spath;
        obj.Info.FileName = sfile;
    end

    [~, ~, ext] = fileparts(sfile);
    if strcmp(ext,'.mat')
        saveAmodel(obj);
    else
        warndlg('Unknown file type')
    end

    if ispath, cd(cpath); end   %return to working folder
end
%%
function saveAmodel(obj)
    %save model setup and results to a mat file as sobj
    spath = obj.Info.PathName;
    sfile = obj.Info.FileName;
    sobj = obj; 
    save([spath,sfile],'sobj');
    clear sobj
end