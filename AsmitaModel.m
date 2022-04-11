classdef AsmitaModel < muiDataSet                       
%
%-------class help---------------------------------------------------------
% NAME
%   AsmitaModel.m
% PURPOSE
%   Class to run Asmita model. Handles run-time initialisation, time
%   stepping and post procressing and holds output in inheritied Data
%   property as a dstable.
% SEE ALSO
%   muiDataSet, ASM_model and ASMinterface for single timestep model
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %first variable in dsproperties list of a new output type
        outType = {'ReachPrism'}; %first reach variable - see PostTimeStep
    end
    
    properties (Transient)
        %Additional properties:     
        Time = 0      %time elapsed from start of run(t=0) in seconds
        DateTime = 0  %time elapsed from Year 0 in seconds
        iStep = 0     %current step number during run
        RunSteps = 0  %number of time steps after checking stability
        delta = 0     %time step in seconds
        outInt = 0    %output interval adjusted for any change in time step
        EleData = []; %array used to store element output at each time step
        RchData = []; %array used to store reach output at each time step
        StepTime      %time to be saved (seconds during run and converted to years)
        MinimumTimeStep = 0.04;  %value used (yrs) as minimum to prompt user
                                 %based on a sp-np cycle (14.6 days)
    end
    
    methods (Access = private)
        function obj = AsmitaModel()                      % << Edit to classname
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = AsmitaModel;                           
            eledsp = elementDSproperties(obj);
            rchdsp = reachDSproperties(obj);
            %check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            %Check additional inputs and notify user what is included
            [message,ok] = AsmitaModel.CheckInput(mobj);
            if ~isempty(message)
                if ~mobj.SupressPrompts  %supress prompts if true
                    button = questdlg(message,'Check','Continue', ...
                        'Quit','Continue');
                    if strcmp(button,'Quit') 
                        delete(obj);
                        return; 
                    end
                end
                %
                if ok<1 
                    warndlg('Not enough components defined to run model');
                    return
                end
            end            

            muicat = mobj.Cases;
%--------------------------------------------------------------------------
% Model code 
%--------------------------------------------------------------------------
            rnpobj = getClassObj(mobj,'Inputs','RunProperties');
            %initialise time step paramters
            tstep = rnpobj.TimeStep; 
            y2s   = mobj.Constants.y2s;    %year to seconds conversion factor
            obj.delta = tstep*y2s;         %time step in seconds
            obj.iStep = 1;                 %initialise model run step
            obj.Time = 0;
            obj.DateTime = rnpobj.StartYear*y2s;
            obj.outInt = rnpobj.OutInterval;
            obj.RunSteps = rnpobj.NumSteps;%may be overwritten in CheckStability
            
            %iniatialise model setup 
            obj = InitialiseModel(obj,mobj,eledsp,rchdsp);
            if isempty(obj), return; end
            %assign the run parameters to the model instance (after
            %InitialiseModel so that ReachID are assigned to Element object
            setRunParam(obj,mobj); %this only saves the core set
            setIncParam(obj,mobj); %add any additional Classes included in run
            
            msg = sprintf('ASMITA processing, please wait');
            hw = waitbar(0,msg);
            %run model
            for jt = 1:obj.RunSteps
                InitTimeStep(obj,mobj,jt)
                RunTimeStep(obj,mobj)
                PostTimeStep(obj,mobj,eledsp,rchdsp);
                waitbar(jt/obj.RunSteps);
                %to report time step during run use the following
                msg = sprintf('ASMITA processing, step %d of %d',...
                                                         jt,obj.RunSteps);
                waitbar(jt/obj.RunSteps,hw,msg);                
            end
            close(hw);

            %now assign results to object properties  
            if rnpobj.StartYear>0
                %code for datetime format
                modeldate = datetime(rnpobj.StartYear,1,1,0,0,0);
                modeltime = modeldate + seconds(obj.StepTime);
                modeltime.Format = eledsp.Row.Format;
            else
                %code for duration format
                modeltime = years(rnpobj.StartYear)+seconds(obj.StepTime);
                eledsp.Row.Format = 'y';
                modeltime.Format = eledsp.Row.Format;
                rchdsp.Row.Format = eledsp.Row.Format;
            end
            %code for calendatDuration format 
            % modeldate = datetime(rnpobj.StartYear,1,1,0,0,0);
            % modeltime = modeldate + seconds(obj.StepTime);
            % modeltime.Format = dsp.Row.Format;
            % modeltime = calendarDuration(datevec(modeltime));
            
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------  
            %assign element data to a dstable
            edst = dstable(obj.EleData{:},'RowNames',modeltime,'DSproperties',eledsp);
            eleobj = getClassObj(mobj,'Inputs','Element');
            elementname = getEleProp(eleobj,'EleName');
            edst.Dimensions.EleName = [elementname;'All'];     %element names used in model run
                       
            %assign metadata about model
            edst.Source = metaclass(obj).Name;
            edst.MetaData = 'Any additional information to be saved';
            
            %assign reach data to a dstable
            rdst = dstable(obj.RchData{:},'RowNames',modeltime,'DSproperties',rchdsp);
            eleobj = getClassObj(mobj,'Inputs','Element');
            idr = unique(getEleProp(eleobj,'ReachID'));
            elementname = getEleProp(eleobj,'EleName');
            rdst.Dimensions.EleName = elementname(idr(idr>0));     %reach names used in model run
                       
            %assign metadata about model
            rdst.Source = metaclass(obj).Name;
            rdst.MetaData = 'Any additional information to be saved';            
            
            %create struct for set of dstables
            dst.Element = edst;  %struct of dstables - fieldnames define Datasets
            dst.Reach = rdst;
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------             
            setDataSetRecord(obj,muicat,dst,'model',[],mobj.SupressPrompts);
            
            %report mass balance on run completion
            asmobj = getClassObj(mobj,'Inputs','ASM_model');
            msg = sprintf('Sediment Balance:  %.3f %%\nWater Balance:  %.3f %%',...
                                        asmobj.SedMbal,asmobj.WatMbal);
            if ~mobj.SupressPrompts  %supress prompts if true
                msgbox(msg,'Run completed')
            end    
        end
%%
        function initialiseModelParameters(mobj)
            %initialise selected parameters that are used in model runs but
            %also for production of graphs and morphological response 
            
            %check whether Reach element types has been modified
            Estuary.setReachTypes(mobj);
            %intitialise transient properties
            Element.initialiseElements(mobj); %used in setEqConcentration
            %initialise equilibrium concentration
            Element.setEqConcentration(mobj);           
            %initial water levels %can be called without initialising
            %AsmitaModel hence set a dummy model obj to initialise time
            obj.Time = 0;
            WaterLevels.setWaterLevels(mobj,obj);
            %initialise dispersion Reach Graph (calls Element.initialiseElements)
            Estuary.initialiseDispersionGraph(mobj); 
            %initialise advection graphs (River and Drift)
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            advobj = Advection.initialiseTransients(mobj);
            if rncobj.IncRiver
                advobj.RiverGraph = Advection.initialiseRiverGraph(mobj);
            end
            %
            if ~isempty(advobj) && rncobj.IncDrift
                advobj.DriftGraph = Advection.initialiseDriftGraph(mobj);
            end
            
            %initialise reaches
            setClassObj(mobj,'Inputs','Reach',[]); %clear any previously created reaches
            Reach.setReach(mobj,true); %sets up reaches and sets properties
            %initialise tidal pumping
            if ~isempty(advobj) && rncobj.IncTidalPumping
                Advection.setTidalPumping(mobj);
            end
            
            setClassObj(mobj,'Inputs','Advection',advobj); %update in case graphs have been added
        end
%%
        function [message,ok] = CheckInput(mobj)
            %Check which inputs are to be included and return message
            %called at run time and as check from main menu
            message = {};  ok=1;
            imes = 1;
            minclasses = mobj.ModelInputs.AsmitaModel;
            msg = 'Core inputs: ';
            for i=1:length(minclasses)
                msg = sprintf('%s %s, ',msg,minclasses{i});
            end
            message{imes} = sprintf('%s are defined\n',msg(1:end-2));  %removes final comma
            
            %check dispersion has been set and that advection balances (if set) 
            estobj = getClassObj(mobj,'Inputs','Estuary');            
            if isempty(estobj.Dispersion)
                imes = imes+1;
                message{imes} = 'Dispersion has not been defined';                
                ok = 0;
            end
            
            %map run conditions to additonal classes required
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            adnclasses = {'River',rncobj.IncRiver;...
                          'Drift',rncobj.IncDrift;...
                          'Saltmarsh',rncobj.IncSaltmarsh;...
                          'Interventions',rncobj.IncInterventions;...
                          'Advection',rncobj.IncTidalPumping;...
                          'CSThydraulics',rncobj.IncDynHydraulics};
                      
            %call to isValidModel has already checked that the minimum
            %set of classes needed to run the model have been initialised
            %now check whether additonal classes to be included are valid
            for j=1:size(adnclasses,1)
                lobj = getClassObj(mobj,'Inputs',adnclasses{j});
                isincluded = adnclasses{j,2};
                if isempty(lobj) && isincluded
                    imes = imes+1;
                    message{imes} = [adnclasses{j,1},' is included but not defined'];
                    ok = 0;
                elseif isincluded
                    imes = imes+1;
                    message{imes} = ['Run includes: ' adnclasses{j,1}];
                end
            end
            
            %check that advection is set and balances (if set) 
            advobj= getClassObj(mobj,'Inputs','Advection'); 
            if (rncobj.IncRiver || rncobj.IncDrift) && isempty(advobj)
                imes = imes+1;
                message{imes} = 'Advection has not been defined';
            elseif rncobj.IncRiver
                okriver = checkMassBalance(advobj,'River');
                if okriver==0
                    imes = imes+1;    
                    message{imes} = 'Mass balance fails for defined River Flows';                
                    ok = 0;
                end            
            end 
        end
    end
%%
    methods (Access = private)
        function obj = InitialiseModel(obj,mobj,eledsp,rchdsp)
             %initialise ASMITA model and internal properties
             
             %NB the interplay between the following functions means that
             %the order in which they are called is important
             %initialise the ASM_model class
             ASM_model.setASM_model(mobj);
             %some initialisation is common to several functions
             AsmitaModel.initialiseModelParameters(mobj);
             
             rncobj = getClassObj(mobj,'Inputs','RunConditions');
             %initialise interventions even if not used (ie set to zero)
             Interventions.initialiseInterventions(mobj);
             %initialise Marsh concentration matrix if marsh included
             if rncobj.IncSaltmarsh
                ok = Saltmarsh.initialiseMarshDepthConc(mobj);
                 if ok<1 %error in concOverMarsh
                     obj = [];
                     return; 
                 end   
             end
             %initialise any river flow or drift equilibrium offset(eqCorV)
             %calls RunConditions.setAdvectionOffset and ASM_model.setDQmatrix
             Element.setEleAdvOffsets(mobj);  %DQmatrix set even if no offset
             %initialise unadjusted equilibirum volumes
             ASM_model.asmitaEqFunctions(mobj);
             %initialise any scaling to initial conditions (kVp)
             Element.setEqScalingCoeffs(mobj);
             %initialise equilibrium conditions
             Element.setEquilibrium(mobj); 
             %constraints and saltmarsh enhanced settling
             ASM_model.setVertExch(mobj,obj);
             %initialise element concentrations
             Element.setEleConcentration(mobj);
             %check time step is small enough
             obj = CheckStability(obj,mobj);
             if isempty(obj)
                 obj = [];
                 return; 
             end
             %initialise mass balance
             ASM_model.MassBalance(mobj,obj);
             
             %write data for initial time step (t=0)
             PostTimeStep(obj,mobj,eledsp,rchdsp); 
        end    
%%
         function InitTimeStep(obj,mobj,jt)
             %initialise model parameters for next time step
             obj.iStep = jt;
             obj.Time = jt*obj.delta;
             obj.DateTime = obj.DateTime+obj.delta;
             tsyear = obj.DateTime/mobj.Constants.y2s;
             %tim = obj.Time/mobj.Constants.y2s;             
             %sprintf('step %g, time %g, year %g',jt,tim,yr)
             rncobj = getClassObj(mobj,'Inputs','RunConditions');
             %update water levels
             WaterLevels.setWaterLevels(mobj,obj);
             %check whether there any forced changes due to interventions
             if rncobj.IncInterventions
                Interventions.setAnnualChange(mobj,obj);
             end
             %update river and drift advection to cater for timeseries flows
             Advection.updateAdvectionGraphs(mobj,tsyear);
             %update reach properties for changes in WL and morphology
             Reach.setReachProps(mobj);
             %assign updated water levels to elements
             Element.setEleWLchange(mobj); 
             %update tidal pumping to reflect morphological changes
             %has an influence if depth, csa or river flow change
             %initialise tidal pumping
             if rncobj.IncTidalPumping
                Advection.setTidalPumping(mobj);
             end
             %update equilibirum volumes due to changes in tidal prism
             Element.setEquilibrium(mobj);
             %update equilibrium concentrations: only needed if cE varied during run
             % Element.setEqConcentration(mobj);
             %set dispersion and advection matrices
             ASM_model.setDQmatrix(mobj,rncobj.Adv2Inc)
             %non-eroding elements + saltmarsh adjustments
             ASM_model.setVertExch(mobj,obj);
             %update element concentrations
             Element.setEleConcentration(mobj);
         end
 %%
        function RunTimeStep(obj,mobj)
            %run model for the time step jt
            ASM_model.asmitaVolumeChange(mobj,obj)  
        end  
%%
        function obj = PostTimeStep(obj,mobj,eledsp,rchdsp)
            %store the results for each time step
            eleobj= getClassObj(mobj,'Inputs','Element'); 

            if obj.iStep==1 || rem(obj.iStep,obj.outInt)==0
                jr = length(obj.StepTime)+1;
                obj.StepTime(jr,1) = obj.Time;
                vname = {eledsp.Variables.Name};  
                M = sign(getEleProp(eleobj,'TransportCoeff'));
                
                for i=1:length(vname)
                    %sorting depends on variable list in DSproperties
                    %element properties
                    eleprop = getEleProp(eleobj,vname{i});
                    if i>3
                        %append mean of all elements
                        meanele = mean(eleprop,1);
                        obj.EleData{i}(jr,:) = [eleprop;meanele];
                    else
                        %append the sum of all elements adjusted for 
                        %water/sediment volumes using sign(n)
                        sumele = sum(M.*eleprop,1);
                        obj.EleData{i}(jr,:) = [eleprop;sumele];
                    end
                end
                
                vname = {rchdsp.Variables.Name};
                for j=1:length(vname)
                    %reach properties
                    rchprop = Reach.getReachProp(mobj,vname{j});
                    obj.RchData{j}(jr,:) = rchprop;  
                end
            end            
        end
%%
        function obj = CheckStability(obj,mobj)
             %Check model stability for given time step
             %If correction for flow equilibrium to be included this needs
             %to be done before checking model stability
             eleobj= getClassObj(mobj,'Inputs','Element'); 
             asmobj = getClassObj(mobj,'Inputs','ASM_model');
             rnpobj = getClassObj(mobj,'Inputs','RunProperties');
             
             vm = getEleProp(eleobj,'MovingVolume');
             %              
             DQ = asmobj.DQ;
             %Estimate concentration in each element
             cnc = getEleProp(eleobj,'EleConcentration');
             check_dt = diag(DQ).*cnc;
             check_dt = min(vm./check_dt)/2;  %the factor of 2 ensures that
                        %dt is substantially less than the stability limit
                       
             obj.RunSteps = rnpobj.NumSteps;                      
             if obj.delta > check_dt
                 y2s = mobj.Constants.y2s;
                 yr_int = round(y2s/check_dt);  %integer number of intervals in a year
                 newtstep = y2s/yr_int;         %adapted time step in seconds
                 dt1 = newtstep/y2s;            %adapted time step in years
                 dt = rnpobj.TimeStep;          %defined time step in years
                 min_dt = obj.MinimumTimeStep; 
                 numtstep = rnpobj.NumSteps*dt/dt1;
                 if dt1<min_dt                     
                     msg = sprintf('Time step check suggests value of %g\nUsing this value will require %g timesteps\nSpecified value is %g',...
                                dt1,numtstep,dt);
                     chktitle = 'Check time step';
                     button = questdlg(msg,chktitle,'Use corrected',...
                         'Use specified','Abort','Abort');
                     switch button
                         case 'Use minimum'
                             dt1 = min_dt;
                         case 'Use specified'
                             dt1 = dt;
                         case 'Abort'
                             obj = [];
                             return;
                     end
                 end
                 obj.outInt = obj.outInt*dt/dt1;    %correct output intervals
                 obj.RunSteps = rnpobj.NumSteps*dt/dt1;%correct number of time steps
                 obj.delta = dt1*y2s;               %revised time step in seconds
             end 
        end
%%
        function setIncParam(obj,mobj)
            %add any additional muiProperty Classes included in run to be
            %saved to RunParam property. 'Element','Estuary','WaterLevels',
            %'RunProperties','RunConditions','EqCoeffParams' classses are 
            %included by default
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncInterventions
                intobj = getClassObj(mobj,'Inputs','Interventions');
                obj.RunParam.Interventions = intobj;
            end
            %
            if rncobj.IncRiver
                rivobj = getClassObj(mobj,'Inputs','River');
                obj.RunParam.River = rivobj;
            end
            %
            if rncobj.IncDrift
                dftobj = getClassObj(mobj,'Inputs','Drift');
                obj.RunParam.Drift = dftobj;
            end
            %
            if rncobj.IncSaltmarsh
                smsobj = getClassObj(mobj,'Inputs','Saltmarsh');
                obj.RunParam.Saltmarsh = smsobj;
            end
            %
            if rncobj.IncDynHydraulics
                cstobj = getClassObj(mobj,'Inputs','CSThydraulics');
                obj.RunParam.CSThydraulics = cstobj;
            end
            %
            if rncobj.IncRiver || rncobj.IncDrift || rncobj.IncDynHydraulics
                advobj = getClassObj(mobj,'Inputs','Advection');
                obj.RunParam.Advection = advobj;                
            end
        end
        
%%
        function dsp = elementDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            %NB variable names in dsp.Variables must match Element 
            % properties for use in PostTimeStep.
            dsp.Variables = struct(...                      
                'Name',{'MovingVolume','FixedVolume','EqVolume',...                           
                           'MovingDepth','FixedDepth','EqDepth',...
                           'BioProdVolume','EleConcentration',...
                           'EleWLchange'},...                           
                'Description',{'Moving Volume','Fixed Volume','Equilibrium Volume',...
                           'Moving Depth','Fixed Depth','Equilibrium Depth',...
                           'Biological Production','Concentration',...                           
                           'Water Level Change'},...
                'Unit',{'m^3','m^3','m^3','m','m','m','m^3','ppm','m'},...
                'Label',{'Volume (m^3)','Volume (m^3)','Volume (m^3)',...
                           'Depth (m)','Depth (m)','Depth (m)',...
                           'Volume (m^3)','Concentration (ppm)',...                           
                           'WL Change (m)'},...
                'QCflag',repmat({'model'},1,9)); 
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'y'},...
                'Label',{'Time (y)'},...
                'Format',{'dd-MMM-yyyy HH:mm:ss'});   %'dd-MMM-yyyy HH:mm:ss'  or 'y'
                                                      %HH for 24 hour clock
            dsp.Dimensions = struct(...    
                'Name',{'EleName'},...
                'Description',{'Element Name'},...
                'Unit',{'-'},...
                'Label',{'Element Name'},...
                'Format',{'-'});  
        end
%%
        function dsp = reachDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            %NB variable names in dsp.Variables must match Element 
            % properties for use in PostTimeStep.
            dsp.Variables = struct(...                      
                'Name',{'ReachPrism','UpstreamPrism',...
                           'UpstreamCSA','RiverFlow','MWlevel',...             
                           'HWlevel','LWlevel','TidalRange'},...                           
                'Description',{'Reach Prism','Tidal Prism',...
                           'Upstream CSA','River Flow','Mean Water Level',...
                           'High Water','Low Water','Tidal Range'},...
                'Unit',{'m^3','m^3','m^2','m/s','mAD','mAD','mAD','m'},...                           
                'Label',{'Volume (m^3)','Volume (m^3)',...
                           'Cross-sectional area (m^2)',...
                           'River flow velocity (m/s)',...
                           'Elevation (mAD)','Elevation (mAD)',...
                           'Elevation (mAD)','Range (m)'},...
                'QCflag',repmat({'model'},1,8)); 
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'y'},...
                'Label',{'Time (y)'},...
                'Format',{'dd-MMM-yyyy HH:mm:ss'});   %'dd-MMM-yyyy HH:mm:ss'  or 'y'
                                                      %HH for 24 hour clock            
            dsp.Dimensions = struct(...    
                'Name',{'EleName'},...
                'Description',{'Element Name'},...
                'Unit',{'-'},...
                'Label',{'Element Name'},...
                'Format',{'-'});  
        end       
    end   
%%
    methods
        function tabPlot(obj,src,mobj) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            tabcb  = @(src,evdat)tabPlot(obj,src,mobj);
            ax = tabfigureplot(obj,src,tabcb,false);            
            
            dst = obj.Data.Element;
            t = dst.RowNames;
            if iscalendarduration(t)
                t = time2num(t);
            end
            vm = dst.MovingVolume;
            vf = dst.FixedVolume;
            ve = dst.EqVolume;
            
            eleobj = getClassObj(mobj,'Inputs','Element');
            n  = getEleProp(eleobj,'TransportCoeff');
            idx = n>0;  %elements that are defined as water volumes
            %total of all water volumens
            vmw = sum(vm(:,idx),2);
            vfw = sum(vf(:,idx),2);
            vew = sum(ve(:,idx),2);
            
            %plot the water elements            
            plot(t,vmw,'DisplayName','Moving')
            hold on
            plot(t,vfw,'DisplayName','Fixed')
            plot(t,vew,'Color',mcolor('green'),'DisplayName','Equilibrium')
            xlabel('Time (y)')
            ylabel('Volume (m^3)')
            title(sprintf('Water volumes for %s',dst.Description))
            hold off
            legend('Location','southeast')
            
            %plot the sediment elements if any
            if any(n<0)
                %totals for all sediment volums
                % vms = sum(vm(:,~idx),2);
                vfs = sum(vf(:,~idx),2);
                ves = sum(ve(:,~idx),2);
                %
                subplot(2,1,1,ax); %convert existing plot to a subplot
                subplot(2,1,2)
                % plot(t,vms,'DisplayName','Moving')
                hold on
                plot(t,vfs,'Color',mcolor('orange'),'DisplayName','Fixed')
                plot(t,ves,'Color',mcolor('green'),'DisplayName','Equilibrium')
                xlabel('Time (y)')
                ylabel('Volume (m^3)')
                title('Sediment volumes')
                hold off
                legend('Location','southeast')
                sgtitle(sprintf('Total change for %s',dst.Description))
            end
        end
    end 
end