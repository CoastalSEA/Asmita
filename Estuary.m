classdef Estuary < muiPropertyUI  
%
%-------class help---------------------------------------------------------
% NAME
%   Estuary.m
% PURPOSE
%   Estuary defines various system properties such as connectivity of
%   elements and reaches, disperison, advection and equilibirum
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
                          'Wind elevation (m)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        Dispersion = []          %dispersin matrix, d (m/s) - uses upper triangle 
        ExternalDisp = []        %array of external exchanges 
    end
    
    properties
        EqRhoCoarse = 0          %coarse fraction equilibrium density (kg/m^3)
        EqRhoFine = 0            %fine fraction equilibrium density (kg/m^3)
        WidthELength = 0         %estuary width e-folding length (m)
        AreaELength = 0          %area convergence length (m)
        WindSpeed = 0            %characteristic annual wind speed (m/s)        
        WindHeight = 10          %wind measurement elevation above msl (m) 
    end    
    
    properties (Dependent)
        EqConcCoarse             %coarse fraction equilibrium concentration(-) 
        EqConcFine               %fine fraction equilibrium concentration(-) 
    end
    
    properties (Transient)       %properties that are time varying
        DispersionGraph               %handle to graph of network (dispersion)
    end

%%   
    methods (Access=protected)
        function obj = Estuary(mobj)     
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
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
            if ~isempty(obj.Dispersion)
                userdata(1,2:end-1) = obj.ExternalDisp(:,1)';
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
                Element.initialiseElements(mobj);                 
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
            
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncDynamicElements
                dispgraph = obj.DispersionGraph;
                [d,exchIn] = graph2matrix(dispgraph);
                dExt = exchIn(:,1);
            else
                dExt = obj.ExternalDisp(:,1);
                d = obj.Dispersion;
            end
            d = d + d';
            D = -d;
            nele = size(d,1);
            for j=1:nele
                D(j,j) = sum(d(j,:))+dExt(j);
            end            
        end
        
%%
        function updateDispersionGraph(mobj)
            %initialise and update dynamic graphs based on modified conditions (if any)            
            obj  = getClassObj(mobj,'Inputs','Estuary');
            if isempty(obj)
                return;
            end
            
            if isempty(obj.DispersionGraph) && ~isempty(obj.Dispersion)
                obj.DispersionGraph = Estuary.initialiseDispersionGraph(mobj);
            else
                %add code to update
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
            if isempty(eleobj)
                warndlg('No elements defined yet');
                return;
            elseif isempty(estobj)
                warndlg('No estuary defined yet');
                return;
            elseif isempty(wlvobj) || abs(wlvobj.SLRrate)==0
                warndlg('SLR rate not set or zero so No Response');
                return;
            end            
            %initialise transient model parameters
            AsmitaModel.intialiseModelParameters(mobj)
            %get required model properties
            y2s = mobj.Constants.y2s;
            eletype = getEleProp(eleobj,'transEleType');
            ich = strcmp(eletype,'Channel');
            finetypes = mobj.GeoType(mobj.FNtypes);
            ifn = matches(eletype,finetypes); %requires v2019b
            vm = getEleProp(eleobj,'MovingVolume');
            sm = getEleProp(eleobj,'MovingSurfaceArea');
            cb = getEleProp(eleobj,'BedConcentration');
            ws = getEleProp(eleobj,'VerticalExchange');
            n  = getEleProp(eleobj,'TransportCoeff');
            cES = estobj.EqConcCoarse;
            cEM = estobj.EqConcFine;
            vhw = Reach.getReachProp(mobj,'HWvolume');
            shw = Reach.getReachProp(mobj,'HWarea');
            dExt = estobj.ExternalDisp;
            %initialise the ASM_model class
            ASM_model.setASM_model(mobj);
            initialiseFlow(rivobj);
            ASM_model.setDQmatrix(mobj,'flow+drift');
            [B,dd] = ASM_model.BddMatrices(mobj);
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
            n_av = mean(abs(n(ifn)));
            tauM = 1/cEM/n_av*(vtot/wso/stot + vtot/delo)/y2s;  %fine
            
            %mask the external elements that are not linked to delta(s)
            extypes = ~contains(mobj.GeoType(mobj.EXtypes),'Delta');            
            iex = matches(eletype,mobj.GeoType(extypes)); %requires v2019b
            slr_max(iex) = 0; slr_max_bio(iex) = 0; mT(iex) = 0;
            
            userData = [slr_max,slr_max_bio,mT];
            Estuary.responseTable(mobj,src,userData,tauS,tauM)
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
            
            if isempty(hx)
                msg1 = 'Maximum rates of sea level rise (mm/year) and morphological response time (Tau)';
                msg2 = 'See Kragtwijk et al (2004) for theoretical basis of estimates';

                helptext = sprintf('%s\n%s',msg1,msg2);
                uicontrol('Parent',src,...
                        'Style','text','String', helptext,...                    
                        'HorizontalAlignment', 'left',...
                        'Units','normalized', 'Position', [0.04 0.85 0.8 0.1],...
                        'Tag','Resptext1');
                    
                msg3 = 'Morphological response times for whole system:';
                helptext = sprintf('%s\n\nChannel/coarse concentration = %.1f years\nFlats/fine concentration = %.1f years',...
                    msg3,tauS,tauM);
                uicontrol('Parent',src,...
                        'Style','text','String', helptext,...                    
                        'HorizontalAlignment', 'left',...
                        'Units','normalized', 'Position', [0.61 0.55 0.4 0.25],...
                        'Tag','Resptext2');
            else
                hxtxt = findobj(hx,'Tag','Resptext2');
                msg3 = 'Morphological response times for whole system:';
                helptext = sprintf('%s\n\nChannel/coarse concentration = %.1f years\nFlats/fine concentration = %.1f years',...
                    msg3,tauS,tauM);
                hxtxt.String = helptext;
            end
        end
    end          

%% ------------------------------------------------------------------------       
% functions that set estuary wide properties
%--------------------------------------------------------------------------
    methods
        function eqConc = get.EqConcCoarse(obj)
            %dependent property derived from EqRhoCoarse
            cn = muiConstants.Evoke;
            eqConc = obj.EqRhoCoarse/cn.SedimentDensity;
        end
%%
        function eqConc = get.EqConcFine(obj)
            %dependent property derived from EqRhoFine
            cn = muiConstants.Evoke;
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
            obj.ExternalDisp(idx) = [];
            obj.Dispersion(idx,:) = [];
            obj.Dispersion(:,idx) = [];
            setClassObj(mobj,'Inputs','Estuary',obj);
        end
    end
end