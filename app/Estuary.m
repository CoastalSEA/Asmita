classdef Estuary < muiPropertyUI  
%
%-------class help---------------------------------------------------------
% NAME
%   Estuary.m
% PURPOSE
%   Estuary defines various system properties such as connectivity of
%   elements and reaches, disperison, advection and equilibrium
%   concentration
% USAGE
%   obj = Estuary.setParamInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Coarse Fraction Equilibirum Density (kg/m^3)',...
                          'Fine Fraction Equilibirum Density (kg/m^3)',...
                          'Width e-folding convergence length (m)',...
                          'Area e-folding convergence length (m)',...
                          'Wind speed (m/s)',...
                          'Wind elevation (m)',...
                          'Flood delta acts as a Reach (1/0)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        Dispersion = []          %dispersion matrix, d (m/s) - uses upper triangle 
        ExternalDisp = []        %array of external exchanges 
        DynamicExchange          %time varying horizontal exchange, used by
                                 %setting exchange in dispersion matrix to NaN
        ExchLinks                %array of vectors defining the from and to
                                 %elements for a dispersion link                    
    end
    
    properties
        EqRhoCoarse = 0          %coarse fraction equilibrium density (kg/m^3)
        EqRhoFine = 0            %fine fraction equilibrium density (kg/m^3)
        WidthELength = 0         %estuary width e-folding length (m)
        AreaELength = 0          %area convergence length (m)
        WindSpeed = 0            %characteristic annual wind speed (m/s)        
        WindHeight = 10          %wind measurement elevation above msl (m) 
        isFDreach = false        %include the Flood Delta as a Reach element
    end    
    
    properties (Dependent)
        EqConcCoarse             %coarse fraction equilibrium concentration(-) 
        EqConcFine               %fine fraction equilibrium concentration(-) 
    end
    
    properties (Transient)       %properties that are time varying
        DispersionGraph          %handle to graph of network (dispersion)
    end

%%   
    methods (Access=protected)
        function obj = Estuary(mobj)     
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
            obj.DynamicExchange = table(0,1e-6,0,'VariableNames',{'Year','Horizontal','Vertical'});
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'Estuary';
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = Estuary(mobj);  
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end  
%%
        function setDispersion(mobj)
            %define or update dispersion properties as an element matrix
            msgtxt = 'Elements not yet defined';
            eleobj = getClassObj(mobj,'Inputs','Element',msgtxt);
            msgtxt = 'Estuary properties not yet defined';
            obj = getClassObj(mobj,'Inputs','Estuary',msgtxt);
            if isempty(eleobj) || isempty(obj), return; end
            
            %assign data for matrix of exchanges
            nele = length(eleobj);
            nelp = nele+2;
            userdata = zeros(nelp,nelp);
            if ~isempty(obj.ExternalDisp) %external exchange has been defined
                userdata(1,2:end-1) = obj.ExternalDisp(:,1)';
            end
            %
            if ~isempty(obj.Dispersion)   %internal exchange has been defined
                userdata(2:end-1,2:end-1) = obj.Dispersion;
            end
            
            %define tableUI for user to edit matrix
            prop = 'Dispersion';
            prompt = 'Enter horizontal exchange between elements (m/s)';
            prompt = sprintf('%s\n(from row element to column element in landward direction) NB: River values are read only.',prompt);
            inoutxt = {'0: Outside';'+: Rivers'};
            [obj.Dispersion,obj.ExternalDisp,obj.DispersionGraph] = ...
                            setmatrix(eleobj,prop,prompt,inoutxt,userdata);
            %assign updated instance
            setClassObj(mobj,'Inputs','Estuary',obj);               
        end     
%%
        function setExchangeRates(mobj)
            %prompt user to define time dependent horizontal exchange rates
            msgtxt = 'Estuary properties not yet defined';
            obj = getClassObj(mobj,'Inputs','Estuary',msgtxt);
            if isempty(obj), return; end
            obj = exchTable(obj);
            setClassObj(mobj,'Inputs','Estuary',obj);
        end
%%
%--------------------------------------------------------------------------
%       INITIALISE
%--------------------------------------------------------------------------
        function [dispersionGraph,nlabel] = initialiseDispersionGraph(mobj)
            %use the Estuary dispersion properties to initialise DispersionGraph
            %when called for a tab plot return graph and labels
            obj  = getClassObj(mobj,'Inputs','Estuary');
            if isempty(obj) || isempty(obj.Dispersion)
                userdata = 0;                
                dispersionGraph = digraph(userdata);
                nlabel = {'No Dispersion Defined'};
            else                 
                Disp = obj.Dispersion;
                ExtDisp = obj.ExternalDisp;
                inoutxt = {'Sea';'Rivers'};
                eleobj  = getClassObj(mobj,'Inputs','Element');
                nodetxt = setnodetext(eleobj,inoutxt);

                dispersionGraph = matrix2graph(Disp,ExtDisp,[],nodetxt);
                nlabel = strcat(num2str(dispersionGraph.Nodes.EleID),...
                                            '-',dispersionGraph.Nodes.Name);
            end
            %when called during run time update instance
            obj.DispersionGraph = dispersionGraph;
            setClassObj(mobj,'Inputs','Estuary',obj);
        end
        
%%
        function setReachTypes(mobj)
            %update the element types that can form reaches if redefined by
            %user
            obj  = getClassObj(mobj,'Inputs','Estuary');
            if obj.isFDreach
                mobj.RCtypes = [1,5];     %Geotypes that define a reach
                mobj.REtypes = [2,3,4,6]; %Geotypes that can belong to a reach 
            else
                mobj.RCtypes = 1;         %Only channel forms a reach (default)
                mobj.REtypes = 2:6;
            end
        end

%%
%--------------------------------------------------------------------------
%       RUNTIME
%--------------------------------------------------------------------------                     
        function [D,dExt] = getDispersion(mobj)
            %get the dispersion properties as an element matrix
            obj  = getClassObj(mobj,'Inputs','Estuary');
            if isempty(obj.DispersionGraph)
                D = 0; dExt = 0;
                return
            end

            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if isempty(rncobj.IncDynamicElements)
                %initialised at run time but not when river or drift flow
                %tab is used before model has been run
                rncobj.IncDynamicElements = false;
            end

            isinitialised = ~any(isnan(obj.DispersionGraph.Edges.Weight));
            if rncobj.IncDynamicElements && isinitialised
                %default is false. set to true if any d=NaN below
                dispgraph = obj.DispersionGraph;
                [d,exchIn] = graph2matrix(dispgraph);
                dExt = exchIn(:,1);
            else
                %use the initial values if the graph has not already
                %been set using updateDispersionGraph
                dExt = obj.ExternalDisp(:,1);
                d = obj.Dispersion;
                idx = find(isnan(d));
                if ~isempty(idx)      
                    rncobj.IncDynamicElements = true; %set dynamic element flag
                    rnpobj = getClassObj(mobj,'Inputs','RunProperties');
                    Estuary.updateDispersionGraph(mobj,rnpobj.StartYear);
                    obj  = getClassObj(mobj,'Inputs','Estuary');
                    
                    dispgraph = obj.DispersionGraph;
                    [d,exchIn] = graph2matrix(dispgraph);
                    dExt = exchIn(:,1);                    
                    [from,to] = ind2sub(size(d),idx);
                    obj.ExchLinks = [from,to];
                end
            end
            d = d + d';
            D = -d;
            nele = size(d,1);
            for j=1:nele
                D(j,j) = sum(d(j,:))+dExt(j);
            end            
        end
        
%%
        function updateDispersionGraph(mobj,tsyear)
            %initialise and update dynamic graphs based on modified conditions (if any)       
            %tsyear is current time-step in years defined in AsmitaModel
            obj  = getClassObj(mobj,'Inputs','Estuary');
            if isempty(obj)
                return;
            end
            
            if isempty(obj.DispersionGraph) && ~isempty(obj.Dispersion)
                %initialise graph if not defined
                obj.DispersionGraph = Estuary.initialiseDispersionGraph(mobj);
            end
            
            if any(isnan(obj.Dispersion),'all')
                %if dispersion matrix includes any NaN values, check
                %whether the dispersion matrix needs to be updated and
                %revise dispersion graph
                Disp = getExchangeRates(obj,tsyear);
                ExtDisp = obj.ExternalDisp;
                inoutxt = {'Sea';'Rivers'};
                eleobj  = getClassObj(mobj,'Inputs','Element');
                nodetxt = setnodetext(eleobj,inoutxt);

                obj.DispersionGraph = matrix2graph(Disp,ExtDisp,[],nodetxt);
            end
            setClassObj(mobj,'Inputs','Estuary',obj);
        end
        
%%
%--------------------------------------------------------------------------
%       OUTPUT
%--------------------------------------------------------------------------   
        function Response(mobj,src,~)
            %present rate of slr and system response times on a tab
            %initialise table for response results
            eleobj = getClassObj(mobj,'Inputs','Element');
            estobj = getClassObj(mobj,'Inputs','Estuary');
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            rivobj = getClassObj(mobj,'Inputs','River');
            rnpobj = getClassObj(mobj,'Inputs','RunProperties'); 
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            cn = getConstantStruct(mobj.Constants);

            if isempty(eleobj) || all([eleobj(:).InitialVolume]==0)
                warndlg('Elements, or element properties, not yet defined');
                return;
            elseif isempty(estobj)
                warndlg('No estuary defined yet');
                return;
            elseif isempty(rnpobj)
                warndlg('No run properties defined yet');
                return;
            elseif isempty(rncobj)
                warndlg('No run conditions defined yet')
                return;
            elseif isempty(wlvobj)
                warndlg('SLR rate not set');
                return;
            else
               startyr = rnpobj.StartYear*cn.y2s;
                newWaterLevels(wlvobj,0,startyr);
                [~,dslr] = getSLR(wlvobj,startyr,startyr); %rate of sea level change (m/y)
                dslr = dslr/cn.y2s;                        %rate of sea level change (m/s)
                if abs(dslr)==0
                    warndlg('SLR is zero so No Response');
                    return;
                end
            end   
            
            %initialise transient model parameters
            AsmitaModel.initialiseModelParameters(mobj)
            %get required model properties
            y2s = mobj.Constants.y2s;
            eletype = getEleProp(eleobj,'transEleType');
            ich = strcmp(eletype,'Channel');
            finetypes = mobj.GeoType(mobj.FNtypes);
            ifn = ismatch(eletype,finetypes); %could replace with matches
            vm = getEleProp(eleobj,'MovingVolume');
            sm = getEleProp(eleobj,'SurfaceArea');
            cb = getEleProp(eleobj,'BedConcentration');
            ws = getEleProp(eleobj,'VerticalExchange');
            n  = getEleProp(eleobj,'TransportCoeff');
            cES = estobj.EqConcCoarse;
            cEM = estobj.EqConcFine;
            vhw = Reach.getReachProp(mobj,'HWvolume');
            shw = Reach.getReachProp(mobj,'HWarea');
            dExt = estobj.ExternalDisp(:,1);
            %initialise the ASM_model class
            ASM_model.setASM_model(mobj);
            if ~isempty(rivobj)
                initialiseFlow(rivobj);
            end
            ASM_model.setDQmatrix(mobj,'flow+drift');
            [B,dd,ok] = ASM_model.BddMatrices(mobj);
            if ok<1, return; end
            %maximum rate of sea level rise for system, mm/year
            slr_max = abs(dd)./cb./sm*y2s*1000;
            %maximum biomass*rate for elements with positive biomass
            Blim = Saltmarsh.limitingBiomass(mobj);
            %maximum rate of sea level rise including any biomass
            slr_max_bio = abs(dd+Blim*vm.*cb)./cb./sm*y2s*1000;
            %system morphological time scale using eigenvalues of the 
            %matrix product B.F
            F = diag(-n./vm);
            BF = B*F;
            lamda = eig(BF);
            mT = -1./lamda/y2s;
            % for system time scale based on Kragtijk et al (2004)
            % sum total volume and plan area at high water
            % should strictly be equilibrium volume for each element
            % use average element properties and horizontal exchange with external
            % environment. The following is a rough approximation when used for a
            % multi-element model
            vtot = sum(vhw); 
            stot = sum (shw);
            wso = mean(ws(ich)); 
            delo = max(dExt);    %this should probably also vary for coarse and fine
            n_av = mean(n(ich));
            tauS = 1/cES/n_av*(vtot/wso/stot + vtot/delo)/y2s; %coarse
            wso = mean(ws(ifn)); 
            n_av = mean(abs(n(ifn)));  %tauM is NaN if ifn=0 no fine types 
            tauM = 1/cEM/n_av*(vtot/wso/stot + vtot/delo)/y2s;  %fine
            
            %mask the external elements that are not linked to delta(s)
            exlist = mobj.GeoType(mobj.EXtypes);
            extypes = ~contains(exlist,'Delta');            
            iex = ismatch(eletype,exlist(extypes)); %could replace with matches
            slr_max(iex) = 0; slr_max_bio(iex) = 0; mT(iex) = 0;
            
            userData = [slr_max,slr_max_bio,mT];
            Estuary.responseTable(mobj,src,userData,tauS,tauM)
        end

%%
        function [Tmin,Tmax] = getMorphTime(mobj)
            %system time scale based on Kragtijk et al (2004) in years
            % sum total volume and plan area at high water
            % should strictly be equilibrium volume for each element
            % use min and max values to estimate range of response time
            eleobj = getClassObj(mobj,'Inputs','Element');
            estobj = getClassObj(mobj,'Inputs','Estuary');
            vhw = Reach.getReachProp(mobj,'HWvolume');
            shw = Reach.getReachProp(mobj,'HWarea');
            ws = getEleProp(eleobj,'VerticalExchange');
            n  = getEleProp(eleobj,'TransportCoeff');
            dExt = estobj.ExternalDisp(:,1);
            cE = estobj.EqConcCoarse;
            y2s = mobj.Constants.y2s;
            vtot = sum(vhw); 
            stot = sum (shw);
            Tmin(1) = 1/cE/max(abs(n))*(vtot/max(ws)/stot + vtot/max(dExt))/y2s; %~minimum
            Tmax(2) = 1/cE/min(abs(n))*(vtot/min(ws)/stot + vtot/min(dExt(dExt>0)))/y2s; %~maximum
        end

%%
        function responseTable(mobj,src,tableData,tauS,tauM)
            %create table for tab to summarise estuary response properties 
            hp = findobj(src,'Type','uipanel');
            hx = findobj(src,'-regexp','Tag','Resptext');
            ht = findobj(src,'Type','uitable');
            delete(ht);
            
            eleobj = getClassObj(mobj,'Inputs','Element');
            rownames = getEleProp(eleobj,'EleName');
            colnames = {'Element Name','slr (no bio)','slr (+bio)','Tau (years)'};
            
            nele = length(rownames);
            nvar = length(colnames);   %number of variables in table
            userdata = cell(nele,nvar);
            for j=1:nele
                userdata{j,1} = rownames{j};
                for k=2:nvar
                    userdata{j,k} = tableData(j,k-1);
                end
            end
            
            if isempty(hp)
                hp = uipanel('Parent',src,...
                    'Units','normalized', ...
                    'Position',[.03 .05 0.57 0.8]);
            end
            uitable('Parent',hp, ...
                'ColumnName', colnames, ...
                'RowName', [], ...
                'ColumnWidth', {100 65 65 65}, ...
                'Data',userdata, ...
                'Units','normalized',...
                'Position',[0,0,1,1]);

            msg0 = sprintf('Morphological response times for whole system (depends on cE, ws, and n, for channel or fine type elements):');
            if isempty(hx)
                msg1 = 'Maximum rates of sea level rise (mm/year) and morphological response time (Tau)';
                msg2 = 'See Kragtwijk et al (2004) for theoretical basis of estimates';

                helptext = sprintf('%s\n%s',msg1,msg2);
                uicontrol('Parent',src,...
                        'Style','text','String', helptext,...                    
                        'HorizontalAlignment', 'left',...
                        'Units','normalized', 'Position', [0.04 0.85 0.8 0.1],...
                        'Tag','Resptext1');
                helptext = sprintf('%s\n\nChannel/coarse concentration = %.1f years\nFlats/fine concentration = %.1f years',...
                    msg0,tauS,tauM);
                uicontrol('Parent',src,...
                        'Style','text','String', helptext,...                    
                        'HorizontalAlignment', 'left',...
                        'Units','normalized', 'Position', [0.61 0.55 0.4 0.25],...
                        'Tag','Resptext2');
            else
                hxtxt = findobj(hx,'Tag','Resptext2');
                helptext = sprintf('%s\n\nChannel/coarse concentration = %.1f years\nFlats/fine concentration = %.1f years',...
                    msg0,tauS,tauM);
                hxtxt.String = helptext;
            end
            
            %Create push button to copy data to clipboard     
            output = num2cell(tableData,1);
            outable = table(output{:},'RowNames',rownames,...
                            'VariableNames',{'slr_nobio','slr_bio','tau'});
            hb = findobj(src,'Tag','uicopy');
            if isempty(hb)
                butpos = [0.75 0.05 0.2 0.05];
                setactionbutton(src,'Copy to clipboard',butpos,...
                               @copydata2clip,'uicopy',...
                              'Copy table content to clipboard',outable);
            else
                hb.UserData = outable;
            end
        end
    end          

%% ------------------------------------------------------------------------       
% functions that set estuary wide properties
%--------------------------------------------------------------------------
    methods
        function eqConc = get.EqConcCoarse(obj)
            %dependent property derived from EqRhoCoarse
            cn = muiConstants.Evoke;  %uses default value of 2650
            eqConc = obj.EqRhoCoarse/cn.SedimentDensity;
        end
%%
        function eqConc = get.EqConcFine(obj)
            %dependent property derived from EqRhoFine
            cn = muiConstants.Evoke;  %uses default value of 2650
            if isempty(obj.EqRhoFine) || obj.EqRhoFine==0
                eqConc = obj.EqConcCoarse;
            else
                eqConc = obj.EqRhoFine/cn.SedimentDensity;
            end
        end
%%       
        function addEleDispersion(obj,mobj)
            %amend dispersion matrix to add an extra element
            % used when element added
            nele = length(obj.ExternalDisp);
            obj.ExternalDisp(nele+1,1) = 0;
            pad = zeros(nele,1);
            obj.Dispersion(nele+1,:) = pad;
            obj.Dispersion(:,nele+1) = [pad;0]';
            setClassObj(mobj,'Inputs','Estuary',obj);
        end       
%%        
        function delEleDispersion(obj,mobj,idx)
            %delete dispersion properties for element(idx)
            % used when element deleted
            Disp = obj.Dispersion;
            ExtDisp = obj.ExternalDisp;
            inoutxt = {'Sea';'Rivers'};
            eleobj  = getClassObj(mobj,'Inputs','Element');
            nodetxt = setnodetext(eleobj,inoutxt);
            dispGraph = matrix2graph(Disp,ExtDisp,[],nodetxt);            
            flowpathID = dispGraph.Nodes.EleID;   %flow path for advection

            idin = find(obj.ExternalDisp(:,1)>0); %element connected to outside environment
            if idx==idin
                inID = find(flowpathID==eleobj(idx).EleID);
                idg = flowpathID(successors(dispGraph,inID)); %next element
                if ~isempty(idg) && idg>0
                    obj.ExternalDisp(idg,1) =  obj.ExternalDisp(idin);
                end
            end

            obj.ExternalDisp(idx,:) = [];
            obj.Dispersion(idx,:) = [];
            obj.Dispersion(:,idx) = [];
            setClassObj(mobj,'Inputs','Estuary',obj);
        end
%%
        function obj = exchTable(obj) 
            %generate UI table for user to edit and add to element interventions
            oldtable = obj.DynamicExchange;
            if isempty(oldtable)
                oldtable = table(0,1e-6,0,'VariableNames',{'Year','Horizontal','Vertical'});
            end
            
            title = 'Dynamic horizontal exchanges';
            header = 'Enter horizontal exchanges for each year required:';
            but.Text = {'Save','Add','Cancel'}; %labels for tab button definition
            newtable = tablefigureUI(title,header,oldtable,true,but,[0.1,0.6]);
            if isempty(newtable), return; end 
            newtable = sortrows(newtable);
            %remove row if year has not been defined
            newtable = rmmissing(newtable,'DataVariables',{'Year'});
            %if user has not edited all NaN values amend to zeros
            newtable = fillmissing(newtable,'constant',0);
            
            obj.DynamicExchange = newtable;
        end
%%
        function d = getExchangeRates(obj,tsyear)
            %dynamic exchange rates are being used - use time step to check
            %whether values need updating.
            d = obj.Dispersion;
            idd = isnan(d);
            idx = find(obj.DynamicExchange.Year<=tsyear,1,'last');            
            d(idd) = obj.DynamicExchange.Horizontal(idx);     
        end
    end
end