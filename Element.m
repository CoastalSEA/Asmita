classdef Element < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   Element.m
% PURPOSE
%   Class for input parameters for some component of the UI application
% USAGE
%   obj = Element.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Initial Volume (m^3)', ...
                          'Initial Surface Area (m^2)',...
                          'Length (m)','Vertical exchange (m/s)',...
                          'Transport coefficient (+/-n)', ...
                          'Bed density (kg/m^3)','Sediment mobility (-)', ...
                          'Tidal Damping (-)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        EleID                   %element id (unique and assigned when instantiated)
                                %NB:does not align with element index once
                                %elements have been deleted
        ReachID = 0             %id of reach numbered in order from mouth 
        EleType                 %element type
        EleName                 %element name
        Erodible = true         %flag to limit erosion of element to initial volume
    end
    
    properties
        InitialVolume = 0       %element volume at start of run (m^3)
        InitialSurfaceArea = 0  %element surface area (m^2)
        Length = 0              %length of element along line of channel (m)
        VerticalExchange = 0    %vertical exchange (m/s)
        TransportCoeff = 3      %transport coefficientn n (-)
        BedDensity = 1650       %bed density(kg/m^3)
        SedMobility = 1         %sediment mobility (used to calibrate tidal pumping)        
        TidalDamping = 1        %option to define tidal amplification/damping factor 
                                %for each element as a/a0. To adjust tidal
                                %range use (2a)*TidalDamping.
    end    

    properties (Transient)
        %Volumes,depths and surface areas during run time
        MovingVolume            %Moving values include water volume change 
        FixedVolume             %Fixed values only include morphological change  
        EqVolume                %Equilbrium volume calculated at run time  
        MovingDepth             %Element depth relative to moving surface        
        FixedDepth              %Element depth relative to fixed surface              
        EqDepth                 %Element equilibrium depth        
        SurfaceArea             %variable surface area (variations within basin area)
        EqSurfaceArea           %Equilibrium surface area calculated at run time
        %other transient element properties
        EqConcentration         %element specific equilbrium concentration
        EleConcentration        %actual concentration in element during run
        transVertExch           %vertical exchange during run including adjustments
        EleWLchange             %water level change that applies to element (hw or lw)
        eqScaling               %scaling of equilibrium relative to initial value
        eqAdvOffSet             %correction to account for advection at t=0
        eqFixedInts             %sum of non-erodible interventions
        BioProdVolume           %cumulative biological production volume
        transEleType            %transient element type (can change during run)
    end    

    properties (Dependent)
        BedConcentration
    end

%%   
    methods (Access=protected)
        function obj = Element(mobj)
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
      function defineElements(mobj)
            %define a set of elements to make up estuary or inlet   
            obj  = getClassObj(mobj,'Inputs','Element');
            if isempty(obj)
                %define numer of elements of each type
                prompt = mobj.GeoType;
                title = 'Number of elements';
                numlines = 1;
                defaultvalues = {'1','1','1','0','0',...
                                 '0','0','0','0','0','0'};                                                 
                %use updated properties to call inpudlg and return new values
                useInp=inputdlg(prompt,title,numlines,defaultvalues);
                if isempty(useInp), return; end %user cancelled

                %assign a type to each element based on user values
                ntype = length(prompt);
                nel_type = str2double(useInp);
                nele = sum(nel_type);
                ele_type = cell(nele,1);
                count = 1;
                for i=1:ntype
                    if nel_type(i)>0
                        for j=1:nel_type(i)
                            etype = char(prompt(i));
                            %rownames{count}=num2str(count);
                            ele_type{count} = etype;
                            count = count+1;
                        end
                    end
                end
                
                %define name and erodibility (and edit Type)
                userdata = table(ele_type,ele_type,true(nele,1));
                userdata = table2cell(userdata); 
                clear obj
                obj = Element(mobj);
                for i=1:nele      
                    obj(i) = Element(mobj); 
                    obj(i).EleID = i;
                    obj(i).EleType = userdata(i,1);
                    obj(i).EleName = userdata(i,2);
                    obj(i).Erodible = userdata{i,3};
                end 
            else
                %elements have been defined so allow user to edit
                obj = mobj.Inputs.Element;
                userdata = [{obj.EleType}',{obj.EleName}',num2cell([obj.Erodible]')];                           
            end
            obj = eleTable(obj,mobj,userdata);
            %check that element names are unique
            while ~isunique({obj(:).EleName})
                hw = warndlg('Element names must be unique','Elements');
                uiwait(hw);
                userdata = [{obj.EleType}',{obj.EleName}',num2cell([obj.Erodible]')];     
                obj = eleTable(obj,mobj,userdata);
            end
            
            setClassObj(mobj,'Inputs','Element',obj);
        end
%%
        function setInput(mobj)
            %gui for user to set Parameter Input values
            msgtxt = 'Use "Define Elements" before defining properties';
            obj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(obj), return; end
            
            eleList = {obj.EleName};
            ok = 1;
            while ok>0
                [idx, ok] = listdlg('Name','Element Properties',...
                'PromptString','Select element','SelectionMode','single',...
                'ListString',eleList);
                if ok==0, continue; end   %Use cancel to quit loop
                %use muiPropertyUI function to generate UI
                lobj = editProperties(obj(idx));
                obj(idx) = lobj;
            end
            setClassObj(mobj,'Inputs','Element',obj);
        end    
%%
        function loadElement(mobj)
            %load element properties from a file
            msgtxt = 'Use "Define Elements" before loading properties';
            obj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(obj), return; end
            
            %reads any tab/space separated file      
            data = readinputdata();
            if isempty(data), return; end     
            
            nele = length(obj);
            ninp = length(data{1});
            if ninp~=nele
                warndlg('Element data in file must match number of elements');
                return;
            end
           
            %values are in order listed in class properties 
            %sort and assign using element ID
            for idx=1:nele                
                idd = find(data{1}==idx);
                obj(idx).InitialVolume = data{2}(idd);
                obj(idx).InitialSurfaceArea = data{3}(idd);
                obj(idx).Length = data{4}(idd);
                obj(idx).VerticalExchange = data{5}(idd);
                obj(idx).TransportCoeff = data{6}(idd);
                obj(idx).BedDensity = data{7}(idd);
                obj(idx).SedMobility = data{8}(idd);
                obj(idx).TidalDamping = data{9}(idd);
            end
            setClassObj(mobj,'Inputs','Element',obj);
        end
%%
        function addElement(mobj)
            %add an element to an existing set of elements
            msgtxt = 'Use "Define Elements" before adding elements';
            obj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(obj), return; end

            eleids = [obj.EleID];
            nelp = length(obj)+1;
            obj(nelp) = Element(mobj);
            obj(nelp).EleID = max(eleids)+1;

            [idx, ok] = listdlg('Name','Element Properties', ...
                'PromptString','Select element', ...
                'SelectionMode','single', ...
                'ListString',mobj.GeoType);
            if ok<1
                warndlg('Cancelled, NO element added')
                return
            end
            
            obj(nelp).EleType = mobj.GeoType{idx};
            obj(nelp).EleName = mobj.GeoType{idx};
            obj(nelp).Erodible = true;
            setClassObj(mobj,'Inputs','Element',obj);
            Element.defineElements(mobj);
             
            %if already created update advection and dispersion
            estobj = getClassObj(mobj,'Inputs','Estuary');
            advobj = getClassObj(mobj,'Inputs','Advection');
            if ~isempty(estobj)
                if ~isempty(estobj.Dispersion)                    
                    addEleDispersion(estobj,mobj);
                end
                if isa(advobj,'Advection')
                    addEleAdvection(advobj,mobj);
                end
            end
            
            %if already created update interventions matrix to match
            intobj = getClassObj(mobj,'Inputs','Interventions');
            if ~isempty(intobj)
                addIntEle(intobj,mobj,obj);
            end
            
            warndlg({'Element ADDED';
              'Do not forget to update Dispersion and Advection'});            
        end              
%%
        function delElement(mobj)
            %delete an element from the model
            msgtxt = 'Use "Define Elements" before deleting elements';
            obj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(obj), return; end
            
            elename = {obj.EleName};
           [idx, ok] = listdlg('Name','Element Properties', ...
                    'PromptString','Select element', ...
                    'SelectionMode','single', ...
                    'ListString',elename);
            if ok<1
                warndlg('Cancelled, NO element deleted')
                return
            end

            %if already created update advection and dispersion
            estobj = getClassObj(mobj,'Inputs','Estuary');
            advobj = getClassObj(mobj,'Inputs','Advection');
            if ~isempty(estobj)
                if ~isempty(estobj.Dispersion)                    
                    delEleDispersion(estobj,mobj,idx);
                end
                %
                if isa(advobj,'Advection')
                    delEleAdvection(advobj,mobj,idx);
                end
            end
            %if already created update interventions matrix to match
            intobj = getClassObj(mobj,'Inputs','Interventions');
            if ~isempty(intobj)
                delIntEle(intobj,mobj,idx);
            end
            warndlg({'Element DELETED';'Properties also removed from:';...
              ' - Dispersion';...
              ' - Advection';...
              ' - Interventions';...
              'Check that the Dispersion and Advection matrices are correctly assigned'});     
            
            %now remove elements (full array used in delEleAdvection to
            %construct graphs and find adjacent elements)
            obj(idx) = [];
            setClassObj(mobj,'Inputs','Element',obj);
        end
%% 
        function initialiseElements(mobj)
            %initialise Transient properties
            msgtxt = 'Unable to initialise because elements not defined';
            obj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(obj), return; end
            
            initdepth = num2cell([obj(:).InitialVolume]./[obj(:).InitialSurfaceArea]); 
            null = num2cell(zeros(1,length(obj)));  
            null2 = num2cell(zeros(2,length(obj)),1);  
            unity = num2cell(ones(1,length(obj)));
            [obj.MovingVolume] = obj(:).InitialVolume;
            [obj.FixedVolume] = obj(:).InitialVolume;
            [obj.EqVolume] = obj(:).InitialVolume;
            [obj.MovingDepth] = initdepth{:};
            [obj.FixedDepth] = initdepth{:};
            [obj.EqDepth] = initdepth{:};
            [obj.SurfaceArea] = obj(:).InitialSurfaceArea;            
            [obj.EqSurfaceArea] = obj(:).InitialSurfaceArea;     
            %initialise default values
            %[obj.EqConcentration]         %set in setEqConcentration
            %[obj.EleConcentration]        %set in setEleConcentration
            [obj.transVertExch] = obj(:).VerticalExchange;     
            [obj.EleWLchange] = null{:};   
            [obj.eqScaling] = unity{:}; 
            %[obj.eqAdvOffSet] = unity{:}; %set in setEleAdvOffsets 
            [obj.eqFixedInts] = null2{:}; 
            [obj.BioProdVolume] =  null{:};
            [obj.transEleType] = obj(:).EleType;                

            setClassObj(mobj,'Inputs','Element',obj);
        end      

%%
        function setEleProp(mobj,varname,prop)
            %assign prop to an element variable name
            msgtxt = 'No elements defined';
            obj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(obj), return; end
            assignum = num2cell(prop);
            [obj.(varname)] = assignum{:};
            setClassObj(mobj,'Inputs','Element',obj);
        end
        
%%                 
        function ok = setEquilibrium(mobj)
            %Set equilibrium volumes for all elements
            %calls AsmitaModel.asmitaEqFunctions(mobj) so that alternatives
            %can be used by overloading in the AsmitaModel class            
            obj = getClassObj(mobj,'Inputs','Element');
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            eqScaling = getEleProp(obj,'eqScaling'); %scaling to initial values
            EqS = getEleProp(obj,'EqSurfaceArea'); 
            %update unadjusted equilibrium volumes to account for changes
            %in tidal prism
            ASM_model.asmitaEqFunctions(mobj);       %updates Element.EqVolume
            %correction for advection            
            eqAdvCor = getEleProp(obj,'eqAdvOffSet');            
            %set adjusted values of equilibrium volume
            for i=1:length(obj)
                EqVol = eqScaling(i)*eqAdvCor(i)*obj(i).EqVolume;
                %make adjustments to volume if intervention are fixed. 
                %Can't do it in ASM_model because EqVolumes are not scaled 
                %to model values. 
                EqVol = EqVol+obj(i).eqFixedInts(1);      %1 - volume change
                if any(EqVol/EqS(i)<0)
                    msgtxt = sprintf('Equilibrium depth in element No.%d is negative, Aborting',i);
                    warndlg(msgtxt); ok = 0; return;
                    % EqVol = 0;
                end                
                [EqVol,~,ok] = checkFlatVolumes(obj,wlvobj,EqVol,0,EqS);
                if ok<1, return; end
                obj(i).EqVolume = EqVol;                
            end
            setClassObj(mobj,'Inputs','Element',obj);
            ok = 1;
        end
        
%%        
        function setEqConcentration(mobj)
            %assign equilibrium concentration to different element types
            obj = getClassObj(mobj,'Inputs','Element');
            %get element types
            nele = length(obj);
            eletype = getEleProp(obj,'transEleType');
            %assign concentration
            finetypes = mobj.GeoType(mobj.FNtypes);
            idx = find(ismatch(eletype,finetypes)); %could be repalced with matches
            
            estobj = getClassObj(mobj,'Inputs','Estuary');
            for i=1:nele
                if any(idx==i)
                    %note: if no fines value defined uses the coarse value
                    %see get.EqConcFine in Estuary class
                    obj(i).EqConcentration = estobj.EqConcFine;
                else
                    obj(i).EqConcentration = estobj.EqConcCoarse;
                end
            end  
            setClassObj(mobj,'Inputs','Element',obj);
        end 
        
%%
        function ok = setEleConcentration(mobj)
            %assign the element concentrations at some instant
            obj = getClassObj(mobj,'Inputs','Element');
            %get element types
            [conc,ok] = ASM_model.asmitaConcentrations(mobj);
            assignum = num2cell(conc);
            [obj.EleConcentration] = assignum{:};
            setClassObj(mobj,'Inputs','Element',obj);
        end
        
%%
        function setEleWLchange(mobj)
            %assign change in high and low water levels at some instant
            obj = getClassObj(mobj,'Inputs','Element');
            %get the Reach element water levels which may have been
            %modified if dynamic hydraulics are being included or the user
            %has specified along channel damping
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncDynHydraulics
                %CSThydraulics computes water levels using the WaterLevels
                %TidalAmp property and a range of river discharges. This 
                %defines a lookup table to obtain the water levels for the
                %dicharge at time, t. CSThydraulics.assignCSTproperties
                %interpolates for each reach element using the dicharge at 
                %time, t, and corrects for any changes in tidal amplitude 
                %and mean sea level at time, t. This is called in 
                %Reach.setHighLowWater which assigns the results for hw, 
                %lw, mwl and tidal range in Reach.setCollectiveProperties.
                %The values at the mouth should be the same as the time
                %step properties in wlvobj.HWaterLevel, etc.
                % water level values derived from CSTmodel
                dHW = Reach.getReachEleProp(mobj,'dHWlevel');
                dLW = Reach.getReachEleProp(mobj,'dLWlevel');
            else
                %apply user defined along channel tidal amplitude damping
                wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
                dHW = wlvobj.dHWchange;  %includes changes due to slr and any
                dLW = wlvobj.dLWchange;  %imposed cycles (ntc etc)
                dMW = wlvobj.dMWchange;  %changes due to slr only
                ampdamp = [obj(:).TidalDamping]; %prescribed amplitude damping
                dHW = (dHW-dMW)*ampdamp+dMW; %apply damping to amplitude only
                dLW = (dLW-dMW)*ampdamp+dMW;
            end

            LWtypes = mobj.GeoType(mobj.LWtypes);   
            %assign the water levels to each element based on LWtypes as
            %defined in Asmita class property.
            nele = length(obj);
            eletype = getEleProp(obj,'transEleType');
            for i=1:nele
                if any(ismatch(LWtypes,eletype(i)))     
                    obj(i).EleWLchange = dLW(i);
                else
                    obj(i).EleWLchange = dHW(i);
                end
            end
            setClassObj(mobj,'Inputs','Element',obj);
        end
%%
        function setEqScalingCoeffs(mobj)
            %function to define scaling for the equilibrium relative to
            %the initial condition          
            obj = getClassObj(mobj,'Inputs','Element');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            
            vm0 = getEleProp(obj,'InitialVolume');            
            ve = getEleProp(obj,'EqVolume');
            if rncobj.ScaleValues
                scale = num2cell(vm0./ve);
                [obj.eqScaling] = scale{:};
            else
                unity = num2cell(ones(1,length(obj)));
                [obj.eqScaling] = unity{:};
            end
            
            %override scaling for saltmarsh elements so that the
            %equilibrium depth is always based on the specified productivity
            eletype = getEleProp(obj,'transEleType');
            ism = contains(eletype,{'Saltmarsh','Storage'});
            if any(ism)
                unity = num2cell(ones(1,sum(ism)));
                [obj(ism).eqScaling] = unity{:};
            end

            setClassObj(mobj,'Inputs','Element',obj);
        end        
%%
        function setEleAdvOffsets(mobj)
            %assign element advection offfset due to river flow or drift
            obj = getClassObj(mobj,'Inputs','Element');          
            nele = length(obj);
            n = getEleProp(obj,'TransportCoeff');           
            %find whether an offset is to be included
            RunConditions.setAdvectionOffset(mobj);
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            %calculate the offset if required    
            
            if strcmp(rncobj.Adv2Offset,'none')
                eqCorV =  ones(nele,1);
            else
                ASM_model.setDQmatrix(mobj,rncobj.Adv2Offset);
                [B,dd] = ASM_model.BddMatrices(mobj);
                eqCorV = (B\dd).^(1./n);                
            end
            %now assign offset to elements 
            assignum = num2cell(eqCorV);
            [obj.eqAdvOffSet] = assignum{:};

            setClassObj(mobj,'Inputs','Element',obj);
            
            %set if 'none' or reset DQ matrix to full value with all advections
            ASM_model.setDQmatrix(mobj,rncobj.Adv2Inc);
        end     
%%
        function setElement(mobj)
            %initialise an empty instance of Element (used in asm_oo2mui)
            obj = Element(mobj);
            setClassObj(mobj,'Inputs','Element',obj);
        end
    end
%%        
        %add other functions to operate on properties as required  
    methods
        function bedConc = get.BedConcentration(obj)
            %get functioon to retrieve the bed concentration as a volume
            %concentration
            cn = muiConstants.Evoke;
            rhow = cn.WaterDensity;
            rhos = cn.SedimentDensity;
            bedConc = (obj.BedDensity-rhow)/(rhos-rhow);
        end

        %%
        function [vm,vf,ok] = checkFlatVolumes(obj,wlvobj,vm,vf,sa)
            %check that tidalflat and saltmarshes are not deeper than tidal
            %range. Should not be needed if coefficients for Ve are 
            %correctly specified and area/volume of element within limit
            fact = 1.0;  %proportion of tidal range for minimum flat            
            TRf = fact*wlvobj.TidalRange;  %factored tidal range
            eletype = getEleProp(obj,'transEleType');
            istype = ismatch(eletype,{'Tidalflat';'Saltmarsh';'DeltaFlat'});         
            mdepths = vm./sa;
            fdepths = vf./sa;
            
            for i=1:length(vm)
                if mdepths(i)>TRf && istype(i)
                    msgtxt = sprintf('Depth in element No.%d is > tidal range, Aborting',i);
                    warndlg(msgtxt); ok = 0; return;
                    % vm(i) = sa(i)*TRf;         %set moving volume to area*tr
                end
                %
                if fdepths(i)>TRf && istype(i)
                    msgtxt = sprintf('Fixed depth in element No.%d > tidal range, Aborting',i);
                    warndlg(msgtxt); ok = 0; return;
                    % vf(i) = sa(i)*TRf;         %set fixed volume to area*tr
                end
            end
            ok = 1;
        end 
%%
        function prop = getEleProp(obj,varname)
            %get a property and return as an element array
            if isempty(obj(1).(varname))
                msg = sprintf('Element property %s not defined',varname);
                warndlg(msg);
                prop = []; return;
            end
            stringtypes = {'EleType';'transEleType';'EleName'};
            if any(strcmp(stringtypes,varname))
                prop = {obj.(varname)}';
            else
                ncol = size(obj(1).(varname),2);
                if ncol>1
                    nrow = size(obj,2); 
                    prop = zeros(nrow,ncol);
                    for i=1:nrow
                        prop(i,:) = obj(i).(varname);
                    end
                else
                    prop = [obj.(varname)]';
                end
            end
        end
%%
        function [eleid,elename,ok] = selectElement(obj)
            %select an element to use for interventions input/modification
            eleList = {obj.EleName};
            elename =[];
            [eleid, ok] = listdlg('Name','Element Change', ...
                'PromptString','Element ID that changes', ...
                'SelectionMode','single', ...
                'ListString',eleList);
            if ok<1, return; end
            elename = eleList{eleid};
        end
%%
        function displayProperties(obj,src)
            %table for summary of element properties.
            %called by muiModelUI.InputTabSummary and overloads super class
            %muiPropertyUI method.            
            rownames = {obj.EleName};
            colnames = {'ID','Element Name','Type','Volume','Surface area','Length', ...
                'ws','n', ...
                'rhob','damp','c(ero)','Erode'};
            props = {'EleID','EleName','EleType','InitialVolume','InitialSurfaceArea',...
                'Length','VerticalExchange','TransportCoeff',...
                'BedDensity','TidalDamping','SedMobility',...
                'Erodible'};
            nele = length(rownames);
            nvar = length(props);   %number of variables in table
            userdata = cell(nele,nvar);
            
            for i=1:nvar
                prop = {obj.(props{i})};
                if isnumeric(prop{1}) && (any([prop{:}]>9999) || ...
                                                     any(abs([prop{:}])<0.001))
                    %when values are e4 or greater change to char e format
                    fmtdata = @(x) sprintf('% .2e',x);
                    prop = cellfun(fmtdata, prop, 'UniformOutput',0);
                end
                for j=1:nele
                    if iscell(prop)
                        userdata{j,i} = prop{j};
                    else
                        userdata{j,i} = prop(j);
                    end
                end
            end
       
            hp = uipanel('Parent',src,...
                'Units','normalized', ...
                'Position',[.01 .01 0.98 0.98]);
            uitable('Parent',hp, ...
                'ColumnName', colnames, ...
                'RowName', {}, ...
                'ColumnWidth', {15,100 55 60 60 60 60 25 45 40 40 40}, ...
                'Data',userdata, ...
                'Units','normalized',...
                'Position',[0,0,1,1]);
        end
    end   
%%    
    methods (Access=private)
        function obj = eleTable(obj,mobj,userdata)
            %populate table of Element Type, Name and erodibility
            typelist = mobj.GeoType;
            data = num2cell(userdata,1);
            oldtable = table(data{:},'VariableNames', {'Type','Name','Erodible'}); 
            oldtable.Properties.UserData.List = {typelist,[],[]};
            title = sprintf('Element Names'); 
            header = 'Edit element Names';  
            but.Text = {'Save','Cancel'}; %labels for tab button definition
            newtable = tablefigureUI(title,header,oldtable,true,but,[0.1,0.6]);
            if isempty(newtable) %user cancelled  
                oldtable.Erodible = [oldtable.Erodible{:}]';
                newtable = oldtable;
            end
            var = newtable{:,1};           
            [obj.EleType] = var{:};            
            var = newtable{:,2};           
            [obj.EleName] = var{:};            
            var = num2cell(newtable{:,3}); 
            [obj.Erodible] = var{:};  
        end     
    end      
end