classdef Advection < matlab.mixin.Copyable
    % Advection class handles different types of advection with functions
    % to create flow graphs and derive advection matrices.
    %
    %----------------------------------------------------------------------
    % AUTHOR
    % Ian Townend
    %
    % COPYRIGHT
    % CoastalSEA, (c) 2016
    %----------------------------------------------------------------------
    %    
    properties               %user defined variables at start of run                      
        RiverFlows = []        %advection matrix, q (m3/s)
        RiverIn = [0,0]            %array of external flows into estuary
        RiverOut = [0,0]           %array of external flows out of estuary
        DriftFlows = []        %matrix of littoral drift flows, qs (m3/s)
        DriftIn = [0,0]            %array of littoral flows into system
        DriftOut = [0,0]           %array of littoral flows out of system
    end
    
    properties (Transient)  %properties that are time varying
        TidalPumpingFlows = [] %tidal pumping matrix, qtp (m3/s)
        TidalPumpingIn = [0,0]     %array of tidal pumping flows into estuary
        TidalPumpingOut = [0,0]    %array of tidal pumping flows out of estuary
        RiverGraph             %handle to graph of river flows (advection)
        DriftGraph             %handle to graph of drift flows (advection)
        QtpGraph               %handle to graph of tidal pumping flows (advection)  
    end
    
    properties (Access=private, Transient)
        InternalAdv = []       %advection matrix, q (m3/s)
        ExternalAdvIn = 0      %matrix of external flows into estuary
        ExternalAdvOut = 0     %matrix of external flows out of estuary        
    end
    
    methods  (Access={?River,?Drift})
        function obj = Advection
        end
    end
    
    methods (Static)
%--------------------------------------------------------------------------
%       SETUP
%--------------------------------------------------------------------------
        function setAdvection(mobj,AdvType)            
            %define or update advection properties as an element matrix
            msgtxt = 'Define Elements before adding Advection';
            eleobj = getClassObj(mobj,'Inputs','Element',msgtxt);
            rivobj = getClassObj(mobj,'Inputs','River');
            dftobj = getClassObj(mobj,'Inputs','Drift');
            %check that elements have already been added and that a river
            %or drift input has been defined
            if isempty(eleobj)
                return;
            elseif isempty(rivobj) &&  strcmp(AdvType,'River')
                warndlg('No river inputs have beend defined')
                return; 
            elseif isempty(dftobj) && strcmp(AdvType,'Drift')
                warndlg('No drift inputs have beend defined')
                return;
            end
            nele = length(eleobj);
            nelp = nele+2;
            
            obj = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj)
                obj = Advection;
            end
            %create or get data and build matrix of exchanges
            userdata = zeros(nelp,nelp);
            obj = setAdvectionType(obj,AdvType);
            if ~isempty(obj.InternalAdv) && size(obj.InternalAdv,1)==nele
                %second condition added to trap incorrectly sized matrix
                %for current number of elements (ie dims have changed)
                userdata(end,2:end-1) = obj.ExternalAdvIn;
                userdata(2:end-1,1) = obj.ExternalAdvOut;
                userdata(2:end-1,2:end-1) = obj.InternalAdv;
            else
                obj.ExternalAdvIn = zeros(nele,1);
                obj.ExternalAdvOut = zeros(nele,1);
            end
            preExtAdvIn = obj.ExternalAdvIn;
            
            %initial labels for table and graph, and setup input matrix
            [qunits,tabletxt,g_flow,flowVar] = assignLabelData(AdvType);  
            
            %get River/Drift inputs - returns struct with 
            % InputEle - id of input element; 
            % okEle - logical flag if correct type of element; 
            % msg - cell array of messages for use in checkSourceInputs;
            % flow - input flows as array the size of River/Drift object array
            src = getSourceProps(g_flow,eleobj,mobj);   
            MatrixInpEle = src.InputEle+1;  %offset due to 'outside'
            userdata(nelp,MatrixInpEle) = src.flow;
                
            %define tableUI for user to edit matrix
            prop = 'Advection';
            prompt = sprintf('Enter advection between elements %s\n(from row element to column element in direction of flow)',qunits);            
            [obj.InternalAdv,AdvIn,AdvOut,obj.RiverGraph] = ...
                            setmatrix(eleobj,prop,prompt,tabletxt,userdata);
            if isempty(AdvIn) || all(AdvIn==0,'all') ||...
                                isempty(AdvOut) || all(AdvOut==0,'all')
                warndlg('One or more exchanges to the outside not defined')
                return; 
            end          
            %setmatrix returns [nx2] exchange arrays. restore to vector
            obj.ExternalAdvIn = AdvIn(:,2);
            obj.ExternalAdvOut = AdvOut(:,1);
            
            %check whether inputs have been changed and update
            checkSourceInputs(obj,mobj,AdvType,flowVar,g_flow,src,preExtAdvIn);
            
            %assign updated flow field to the specified flow type
            obj = addAdvectionProps(obj,AdvType);   
            setClassObj(mobj,'Inputs','Advection',obj);
            
            %--------------------------------------------------------------
            %nested function to initialise labels for table and graph
            %--------------------------------------------------------------
            function [qunits,tabletxt,g_flow,flowVar] = ...
                                                 assignLabelData(AdvType)
                %initialise labels for table and graph                
                switch AdvType
                    case 'River'
                        if ~isempty(rivobj)
                            g_flow = rivobj;
                            flowVar = 'RiverFlow';    
                        end
                        qunits = '(m3/s)';
                        tabletxt = {'0: Sea';'+: Rivers'};
                    case 'Drift'
                        if ~isempty(dftobj)
                            g_flow = dftobj;
                            flowVar = 'DriftRate';
                            
                        end
                        qunits = '(m3/year)';
                        tabletxt = {'>: Sink';'>: Source'};
                    otherwise
                        qunits = ''; tabletxt = {'',''};
                        g_flow = []; flowVar = '';
                end
            end
            %--------------------------------------------------------------
            %nested function to check for changes in source inputs
            %--------------------------------------------------------------
            function checkSourceInputs(obj,mobj,AdvType,flowVar,...
                                                   g_flow,src,preExtAdvIn)                 
                %check whether inputs have been changed and update
                idx = diff([preExtAdvIn,obj.ExternalAdvIn],1,2);
                idinputs = find(idx~=0);
                if ~isempty(src.InputEle)        %if some source already exists
                    idnewinput = setdiff(idinputs,src.InputEle);
                    if ~isempty(idnewinput)
                        for i=1:length(idnewinput)   %add new sources                        
                            %check that source is assigned to an element that exists   
                            ok = src.okEle(idnewinput(i));
                            if ok==0 %catch changes to elements that are not valid
                                if obj.ExternalAdvIn(idnewinput(i))>0 %only report when a non-zero value has been input (ie not corrections back to zero)
                                    msgtxt = sprintf(src.msg{1},idnewinput(i));
                                    warndlg(msgtxt);
                                end
                            else
                                %add a new source
                                addSource(obj,mobj,AdvType,idnewinput(i));
                                msgtxt = sprintf(src.msg{2},idnewinput(i));
                                warndlg(msgtxt);
                            end
                        end
                    end
                    
                    for i=1:length(src.InputEle) 
                        %update existing sources
                        g_flow(i).(flowVar) =  obj.ExternalAdvIn(src.InputEle(i));
                    end

                else                         %if no input source exist
                     if ~isempty(idinputs)
                         warndlg('Inputs have been added and need to be defined in Setup');
                         for i=1:length(idinputs)
                             addSource(obj,mobj,AdvType,idinputs(i));
                         end
                     end
                end
            end
        end
        
%%       
%--------------------------------------------------------------------------
%       INITIALISE
%--------------------------------------------------------------------------
        function [flowGraph,nlabel] = initialiseRiverGraph(mobj)
            %use initial River advection properties to initialise RiverGraph
            obj  = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj) || isempty(obj.RiverFlows) || all(obj.RiverIn(:)==0)
                userdata = 0;
                flowGraph = digraph(userdata);
                nlabel = {'No River Flows Defined'};
            else
                Element.initialiseElements(mobj);
                Adv = obj.RiverFlows;
                AdvIn = obj.RiverIn;
                AdvOut = obj.RiverOut;
                inoutxt = {'Sea';'River(s)'};
                eleobj  = getClassObj(mobj,'Inputs','Element');
                nodetxt = setnodetext(eleobj,inoutxt);
                flowGraph = matrix2graph(Adv,AdvIn,AdvOut,nodetxt);
                nlabel = strcat(num2str(flowGraph.Nodes.EleID),...
                                                '-',flowGraph.Nodes.Name);
                                            
                %find if river inputs have changed and update if required
                rivobj = getClassObj(mobj,'Inputs','River');                                        
                j = zeros(length(rivobj),1);                 
                for i=1:length(rivobj) %find souces defined by a timeseries
                   j(i) = ~isempty(rivobj(i).RiverTSC); 
                end               
                %
                if sum(j)>0 %update if timeseries input
                    obj.RiverGraph = flowGraph; %update the model instance
                    setClassObj(mobj,'Inputs','Advection',obj);
                    %
                    rnpobj = getClassObj(mobj,'Inputs','RunProperties');
                    Advection.updateAdvectionGraphs(mobj,rnpobj.StartYear);
                    flowGraph = obj.RiverGraph; %reassign in case it has been updated
                else
                    initialiseFlow(rivobj);     %initialise the transient values
                end                                         
                
                %check mass balance of advective flows
                ok = checkMassBalance(obj,'River');
                if ok<1
                    warndlg('Mass balance fails for defined Advection');
                end                
            end
        end
      
%%
        function [driftGraph,nlabel] = initialiseDriftGraph(mobj)
            %use initial Drift advection properties to initialise DriftGraph
            obj  = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj) || isempty(obj.DriftFlows) || all(obj.DriftIn(:)==0)
                userdata = 0;
                driftGraph = digraph(userdata);
                nlabel = {'No Drift Flows Defined'};
            else
                Element.initialiseElements(mobj);
                Qs = obj.DriftFlows;
                QsIn = obj.DriftIn;
                QsOut = obj.DriftOut;
                inoutxt = {'Downdrift','Updrift'};
                eleobj  = getClassObj(mobj,'Inputs','Element');             
                nodetxt = setnodetext(eleobj,inoutxt);
                driftGraph = matrix2graph(Qs,QsIn,QsOut,nodetxt);
                nlabel = strcat(num2str(driftGraph.Nodes.EleID),...
                    '-',driftGraph.Nodes.Name);
                
                %find if drift inputs have changed and update if required
                dftobj = getClassObj(mobj,'Inputs','Drift');                                    
                j = zeros(length(dftobj),1); 
                for i=1:length(dftobj) %find souces defined by a timeseries
                   j(i) = ~isempty(dftobj(i).DriftTSC); 
                end    
                %
                if sum(j)>0  %update if timeseries input
                    obj.DriftGraph = driftGraph; %update the model instance
                    setClassObj(mobj,'Inputs','Advection',obj);
                    rnpobj = getClassObj(mobj,'Inputs','RunProperties');
                    Advection.updateAdvectionGraphs(mobj,rnpobj.StartYear);
                    driftGraph = obj.DriftGraph;
                else
                    initialiseFlow(dftobj);
                end                                           
            end           
        end
        
%%
        function [qtpGraph,nlabel] = initialiseQtpGraph(mobj)
            %use the tidal pumping properties to initialise qtpGraph
            userdata = 0;
            qtpGraph = digraph(userdata);
            obj  = getClassObj(mobj,'Inputs','Advection');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if isempty(rncobj)
                nlabel = {'Need to define Run Parameters>Conditions'};
            elseif ~rncobj.IncTidalPumping         
                nlabel = {'Include Tidal Pumping condition not set'};
            elseif isempty(obj.RiverFlows)
                nlabel = {'No river flows defined so there is no tidal pumping'};                
            else
                %intitialise transient properties - calls Advection.setTidalPumping
                AsmitaModel.initialiseModelParameters(mobj);                 
                if isempty(obj.TidalPumpingFlows) %if still not defined
                    msg = sprintf('Check parameters for tidal pumping have been defined');
                    nlabel = {msg};
                    return;
                end
                Qtp = obj.TidalPumpingFlows;
                QtpIn = obj.TidalPumpingIn;
                QtpOut = obj.TidalPumpingOut;
                inoutxt = {'Sea';'River(s)'};
                eleobj  = getClassObj(mobj,'Inputs','Element');
                nodetxt = setnodetext(eleobj,inoutxt);
                qtpGraph = matrix2graph(Qtp,QtpIn,QtpOut,nodetxt);
                nlabel = strcat(num2str(qtpGraph.Nodes.EleID),...
                                                '-',qtpGraph.Nodes.Name);
            end
        end
       
%%
%--------------------------------------------------------------------------
%       RUNTIME
%--------------------------------------------------------------------------
        function obj = initialiseTransients(mobj)
            %initialise transient Advection properties at start of model run
            obj = getClassObj(mobj,'Inputs','Advection');
            obj.TidalPumpingFlows = []; %tidal pumping matrix, qtp (m3/s)
            obj.TidalPumpingIn = 0;     %array of tidal pumping flows into estuary
            obj.TidalPumpingOut = 0;    %array of tidal pumping flows out of estuary
            obj.RiverGraph = [];        %handle to graph of river flows (advection)
            obj.DriftGraph = [];        %handle to graph of drift flows (advection)
            obj.QtpGraph = [];          %handle to graph of tidal pumping flows (advection)  

            obj.InternalAdv = [];       %advection matrix, q (m3/s)
            obj.ExternalAdvIn = 0;      %matrix of external flows into estuary
            obj.ExternalAdvOut = 0;     %matrix of external flows out of estuary             
        end
%%
        function [Q,qIn,qOut] = getAdvectionFlow(mobj,AdvType)
            %get the advection properties as an element matrix
            AdvGraphName = strcat(AdvType,'Graph');     
            obj  = getClassObj(mobj,'Inputs','Advection');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if isempty(obj) || isempty(obj.(AdvGraphName)) || ...
                               isempty(obj.(AdvGraphName).Nodes)
                Q = 0; qIn = 0; qOut = 0;
                return;
            end
            
            if isempty(rncobj.IncDynamicElements)
                %initialised at run time but not when river or drift flow
                %tab is used before model has been run
                rncobj.IncDynamicElements = false;
            end

            %check whether river or drift use timeseries
            isdynamic = false;
            switch AdvType
                case 'Drift'
                    if rncobj.IncDrift
                        %only include time series if drift is to be included
                        isdynamic = rncobj.IncDriftTS;
                    end
                case 'River'
                    if rncobj.IncRiver
                        %only include time series if river is to be included
                        isdynamic = rncobj.IncRiverTS;
                    end    
            end
            
            %get flow properties using graph if dynamic or properties if not
            nele = length(getClassObj(mobj,'Inputs','Element'));
            if rncobj.IncDynamicElements || isdynamic
                AdvGraph = obj.(AdvGraphName);
                [q,exchIn,exchOut] = graph2matrix(AdvGraph,nele);
                qIn = exchIn(:,2);
                qOut = exchOut(:,1);
            else 
                obj = setAdvectionType(obj,AdvType);
                qIn = obj.ExternalAdvIn;
                qOut = obj.ExternalAdvOut;
                q = obj.InternalAdv;
            end
            
            if strcmp(AdvType,'Drift')
                %convert drift rates from m3/yr to concentration equivalent
                %flow rate
                [q,qIn,qOut] = getDriftFlow(obj,mobj,q,qIn,qOut);
                if isnan(q),  Q = NaN;  return; end %
            end
            
            %set up flow matrix
            Q = advMatrix(obj,q,qOut,nele);
        end 
   
%% 
        function updateAdvectionGraphs(mobj,tsyear)
            %initialise and update dynamic graphs based on modified conditions (if any)            
            obj = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj)
                return;
            end
            
            rivobj  = getClassObj(mobj,'Inputs','River');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            %only check need to update if river is included
            if rncobj.IncRiver && ~isempty(obj.RiverGraph) && ~isempty(rivobj)
                obj.RiverGraph = getflowpath(obj,mobj,rivobj,tsyear);                
            end
            
            dftobj  = getClassObj(mobj,'Inputs','Drift');
            %only check need to update if drift is included
            if rncobj.IncDrift && ~isempty(obj.DriftGraph) && ~isempty(dftobj)
                obj.DriftGraph = getflowpath(obj,mobj,dftobj,tsyear);
            end

            setClassObj(mobj,'Inputs','Advection',obj);
             
            %--------------------------------------------------------------
            %nested function to construct the relevant flowgraph
            %--------------------------------------------------------------
            function g_flowpath = getflowpath(obj,mobj,srcobj,tsyear)
                %get relevant flowgraph (river or drift based on srcobj)
                if isa(srcobj,'River')
                    GraphName = 'RiverGraph';
                    TSname = 'RiverTSC';
                elseif isa(srcobj,'Drift')
                    GraphName = 'DriftGraph';
                    TSname = 'DriftTSC';
                else
                    return;
                end
                g_flowpath = obj.(GraphName); %existing definition
                j = zeros(length(srcobj),1); 
                for is=1:length(srcobj) %find sources defined by a timeseries
                   j(is) = ~isempty(srcobj(is).(TSname)); 
                end                
                if sum(j)>0  %timeseries present so update flowgraph
                    for is=1:length(j) %loop to consider each source                        
                        if j(is)>0
                            %interpolate time series to get new input value
                            if isa(srcobj,'Drift') 
                                [inflow,inEleID] = Drift.getDriftTSprop(mobj,...
                                                            is,tsyear);                                                
                            else
                                [inflow,inEleID] = River.getRiverTSprop(mobj,...
                                                            is,tsyear);
                            end
                            %
                            if length(inflow)>1 %drifts vary along path %   CHECK IF CONDITION
                                g_flowpath = Advection.newDriftPath(g_flowpath,...
                                                        inflow,inEleID);
                            else %source value dictates flow along path
                                inEleGraph = find(g_flowpath.Nodes.EleID==inEleID)-1;
                                isbalance = true;   %ie update whole network 
                                %get the matrix and exchages for the input graph
                                [~,InFlow] = graph2matrix(g_flowpath);
                                if any(InFlow(inEleGraph,2)~=inflow) 
                                    %update if any input has changed
                                    InFlow(inEleGraph,2) = inflow; %assigne changes to [nx2] input array
                                    g_flowpath = rescale_graph(g_flowpath,...
                                                        InFlow,isbalance);
                                end
                            end
                        end
                    end
                end 
            end
        end
%%
        function g_flowpath = newDriftPath(g_flowpath,newDrift,inEleID)
            %reassign new flows to flowgraph when there are multiple sources. 
            %newDrift - [nx1] vector of drifts to be updated
            %inEleID - [nx2] array of start and end node for each drift input
            %          [0,j] = input and [j,0] = output
            %mass balance is not checked because arbitrary inputs can be
            %specified for any edge.
            flowpathID = g_flowpath.Nodes.EleID;
            newEdgeWeight = g_flowpath.Edges.Weight;
            %find Edge for each dynamic drift input and assign update
            for j=1:length(newDrift)
                if inEleID(j,1)==0      %an input to node (j,2)                      
                    outID = find(flowpathID==inEleID(j,2));
                    upstreamids = predecessors(g_flowpath,outID);
                    inID = upstreamids((flowpathID(upstreamids)==0));
                elseif inEleID(j,2)==0 %an output from node (j,1)
                    inID = find(flowpathID==inEleID(j,1));
                    dnstreamids = successors(g_flowpath,inID);
                    outID = dnstreamids((flowpathID(dnstreamids)==0));
                else                   %exchange between (j,1) and (j,2)
                    inID = find(flowpathID==inEleID(j,1));
                    outID = find(flowpathID==inEleID(j,2));
                end
                ide = findedge(g_flowpath,inID,outID);
                newEdgeWeight(ide) = newDrift(j);
            end
            g_flowpath.Edges.Weight = newEdgeWeight;
        end
%%
        function setTidalPumping(mobj)
            %assign the tidal pumping discharge to the flow network
            % Calculate the advective flux up estuary due to tidal pumping
            % Based on modified version of equation (18) proposed by Friedrichs,
            % Armbrust and de Swart, 1998.
            obj  = getClassObj(mobj,'Inputs','Advection');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            estobj = getClassObj(mobj,'Inputs','Estuary');
            
            if ~rncobj.IncTidalPumping || ~isa(obj,'Advection')                
                return;
            elseif isempty(wlvobj) || isempty(estobj)
                return;
            else
                if isempty(obj.RiverGraph)     %intitialise RiverGraph
                    obj.RiverGraph = Advection.initialiseRiverGraph(mobj);
                end                
            end

            %get reach based tidal pumping discharge
            [qtp,qtp0] = Advection.getTidalPumpingDischarge(mobj);

            if any(qtp>0)
                g_landward = tidalpumpingGraph(obj,mobj,qtp,qtp0);
                %assign instance values to internal properties
                nele = size(obj.RiverFlows,1);
                [qtpFlows,qtpIn,qtpOut,~] = graph2matrix(g_landward,nele);
                obj.TidalPumpingIn = qtpIn;
                obj.TidalPumpingOut = qtpOut;
                obj.TidalPumpingFlows = qtpFlows;
                obj.QtpGraph = g_landward;
                %assign updated flow field to the tidal pumping properties
                setClassObj(mobj,'Inputs','Advection',obj);
            end
        end  
%%
        function [qtp,qtp0] = getTidalPumpingDischarge(mobj)
            %compute tidal pumping for each reach with a river flow based on tide and river properties
            wlvobj  = getClassObj(mobj,'Inputs','WaterLevels');
            aM4  = wlvobj.M4amplitude; %amplitude of M4 tidal constituent (m)        
            phM4 = wlvobj.M4phase;     %phase of M4 tidal constituent (m)            
            zHW0 = wlvobj.HWaterLevel; %High Water Level at mouth at time t         
            zLW0 = wlvobj.LWaterLevel; %Low Water Level at mouth at time t
            Tp   = wlvobj.TidalPeriod; %tidal period (hrs)            

            eleobj = getClassObj(mobj,'Inputs','Element');
            eLe = getEleProp(eleobj,'Length');
            ws = getEleProp(eleobj,'transVertExch');
            cero = getEleProp(eleobj,'SedMobility');
            
            estobj = getClassObj(mobj,'Inputs','Estuary');
            eLw  = estobj.WidthELength;%estuary width e-folding length (m)
            eLa  = estobj.AreaELength; %estuary area e-folding length (m)
            
            idrch = Reach.getReachProp(mobj,'ReachChannelID');
            Hhw = Reach.getReachProp(mobj,'HWvolume')./... %depth at high water (m)
                           Reach.getReachProp(mobj,'HWarea');
            Hlw = Reach.getReachProp(mobj,'LWvolume')./... %depth at low water (m)
                           Reach.getReachProp(mobj,'LWarea');

            %derived variables needed to calculate qtp              
            aef = (Hhw-Hlw)/2;        %effective tidal range (m)
            aef(aef<=0) = 0.0001;     %prevent aef from going negative for channels with very large tidal flats 
            aef = mean(aef);          %estuary mean value effective tidal range(m)
            Hav = (Hhw+Hlw)/2;        %average depth of elements(m)
            Hav = mean(Hav);          %average depth of estuary (channel only)(m)
            sigma = 2*pi/Tp/3600;     %tidal frequency (1/s)
            amp = (zHW0-zLW0)/2;      %tidal amplitude at mouth (m)
            M4M = aM4/amp;            %ratio of M4 to tidal amplitude at mouth
            %
            U=(sigma*amp.*eLw)./Hav;         
            cd=(0.4./(log(Hav./0.0003)-1)).^2;
            F=(8/(3*pi).*cd.*U./Hav);
            kno=F.*U./(9.81.*amp);    %tidal wave number (1/m)
            kno(isnan(kno)) = 0;
            
            %cumulative length from mouth, xi, to upstream end of element
            xi = Reach.getReachProp(mobj,'CumulativeLength');
            %area at upstream section based on average of adjacent elements
            upCSA = Reach.getReachProp(mobj,'UpstreamCSA'); %best validation in Severn case
            dnCSA = upCSA.*exp(eLe(idrch)/eLa);             %but sensitive to eLa and eLw in U
            %alternative basis for area - smoother description of CSA but
            %delays increase in conc to much further upstream in Severn case
            % CSA = Reach.getReachProp(mobj,'ReachCSA');
            % upCSA = CSA.*exp(-eLe(idrch)/2/eLa);
            % dnCSA = CSA.*exp(eLe(idrch)/2/eLa);

            %derive river flow speed from river discharge and local(upstream)csa
            uriv = Reach.getReachProp(mobj,'RiverFlow');
            %
            % tpump = (3/2.cero.U.(aeff.k.x/H-(a2/a).sin(phM4)+1/4.Tc.U/Lw)-uriv*(3/2.cero-1)).Axsec  [eq.1]
            % U  = sig*Lw*a/H; Tc = H/ws;
            % cero is defined as alpha' in the source paper and manual
            % compute tidal pumping flux per unit area. Re-write eqn.1 per unit area as:
            % tpump'=1.5*cero*sig*eLw*amp/Hav*(aef*kno*xi/Hav-M4M*sin(phM4)+0.25*sig*amp/ws)-uriv*(3/2.cero-1)
            term1 = 3/2*cero(idrch).*U;
            term2 = aef.*kno.*xi./Hav;
            term3 = M4M*sin(phM4*pi/180); 
            term4 = 0.25*sigma*amp./ws(idrch);
            term5 = uriv.*(3/2*cero(idrch)-1);
            %tidal pumping at the upstream end of each reach 
            qtp = floor((term1.*(term2-term3+term4)-term5).*upCSA);
            %tidal pumping at the downstream end of each reach  
            term2 = aef.*kno.*(xi-eLe(idrch))./Hav;
            term5 = [term5(1);term5(1:end-1)]; %offset to get downstream velocity      
            qtp0 = floor((term1.*(term2-term3+term4)-term5).*dnCSA); %only mouth values used in tidalpumpingGraph
        end

%%       
        function inputSummary(mobj,src,~)
            %setup table of flow inputs for Flows>Input Summary tab
            ht = findobj(src,'Type','uitable');
            delete(ht);
            ht = findobj(src,'Tag','IStext');
            delete(ht);
            %Headers
            uicontrol('Parent',src,...
                    'Style','text','String','River properties',...                    
                    'HorizontalAlignment', 'left',...
                    'Units','normalized', 'Position', [0.05 0.9 0.4 0.05],...
                    'Tag','IStext'); 
            uicontrol('Parent',src,...
                    'Style','text','String', 'Drift properties',...                    
                    'HorizontalAlignment', 'left',...
                    'Units','normalized', 'Position', [0.05 0.45 0.4 0.05],...
                    'Tag','IStext');
            %get tables for River and Drift  
            rivobj  = getClassObj(mobj,'Inputs','River');
            if ~isempty(rivobj)
                riverPropertiesTable(rivobj,src,mobj);
            end
            dftobj  = getClassObj(mobj,'Inputs','Drift');
            if ~isempty(dftobj)
                driftPropertiesTable(dftobj,src,mobj);
            end
        end       
    end
    
%%
    methods
        function addEleAdvection(obj,mobj)
            %amend advection matrix to add an extra element
            if ~isempty(obj.RiverFlows)
                obj = setAdvectionType(obj,'River');
                obj = updateAdvectionArrays(obj);
                obj = addAdvectionProps(obj,'River');
            end
            %
            if ~isempty(obj.DriftFlows)
                obj = setAdvectionType(obj,'Drift');
                obj = updateAdvectionArrays(obj);
                obj = addAdvectionProps(obj,'Drift');
            end
            %
            setClassObj(mobj,'Inputs','Advection',obj);        
            %--------------------------------------------------------------
            function obj = updateAdvectionArrays(obj)
                nele = length(obj.ExternalAdvIn);
                obj.ExternalAdvIn(nele+1,1) = 0;
                obj.ExternalAdvOut(nele+1,1) = 0;
                pad = zeros(nele,1);
                obj.InternalAdv(nele+1,:) = pad;
                obj.InternalAdv(:,nele+1) = [pad;0]';
            end
        end
        
%%        
        function delEleAdvection(obj,mobj,idx)
            %delete advection properties for element(idx)
            if ~isempty(obj.RiverFlows)
                Adv = obj.RiverFlows;
                AdvIn = obj.RiverIn;
                AdvOut = obj.RiverOut;
                inoutxt = {'Sea';'River(s)'};
                eleobj  = getClassObj(mobj,'Inputs','Element');
                nodetxt = setnodetext(eleobj,inoutxt);
                flowGraph = matrix2graph(Adv,AdvIn,AdvOut,nodetxt);
                obj = setAdvectionType(obj,'River');
                obj = deleteAnAdvection(obj,eleobj,flowGraph,idx);
                obj = setAdvectionProps(obj,'River',idx);
            end
            %
            if ~isempty(obj.DriftFlows)
                Qs = obj.DriftFlows;
                QsIn = obj.DriftIn;
                QsOut = obj.DriftOut;
                inoutxt = {'Downdrift','Updrift'};
                eleobj  = getClassObj(mobj,'Inputs','Element');             
                nodetxt = setnodetext(eleobj,inoutxt);
                driftGraph = matrix2graph(Qs,QsIn,QsOut,nodetxt);
                obj = setAdvectionType(obj,'Drift');
                obj = deleteAnAdvection(obj,eleobj,driftGraph,idx);
                obj = setAdvectionProps(obj,'Drift',idx);
            end
            %
            %assign updated flow field to the specified flow type
            setClassObj(mobj,'Inputs','Advection',obj);

            %--------------------------------------------------------------
            function obj = deleteAnAdvection(obj,eleobj,flowGraph,idx)
                %check if end member is being deleted and move input/output
                %to adjacent element
                flowpathID = flowGraph.Nodes.EleID; %flow path for advection
                idin = find(obj.ExternalAdvIn>0); %element connected to outside environment
                if idx==idin
                    inID = find(flowpathID==eleobj(idx).EleID); %incoming path
                    idg = flowpathID(successors(flowGraph,inID));
                    if ~isempty(idg) && idg>0
                        obj.ExternalAdvIn(idg) =  obj.ExternalAdvIn(idin);
                    end
                end
                obj.ExternalAdvIn(idx) = [];

                %nele = length(obj.ExternalAdvOut);
                idout = find(obj.ExternalAdvOut>0); %element connected to outside environment
                if idx==idout
                    outID = find(flowpathID==eleobj(idx).EleID); %outgoing path
                    idg = flowpathID(predecessors(flowGraph,outID));
                    if ~isempty(idg) && idg>0
                        obj.ExternalAdvOut(idg) =  obj.ExternalAdvOut(idout);
                    end
                end
                obj.ExternalAdvOut(idx) = [];

                obj.InternalAdv(idx,:) = [];
                obj.InternalAdv(:,idx) = [];
            end
        end

%%
        function [qs,qsIn,qsOut] = getDriftFlow(obj,mobj,QS,QSIn,QSOut)            
            %adjust littoral drift rates to equivalent flow rates
            eleobj = getClassObj(mobj,'Inputs','Element');
            estobj = getClassObj(mobj,'Inputs','Estuary');
            cE = getEleProp(eleobj,'EqConcentration');
            nele = length(eleobj);
            kCeI = River.getRiverProp(mobj,'tsRiverConc')./cE; 
            y2s = mobj.Constants.y2s;
            cES = estobj.EqConcCoarse;
            [D,dExt] = Estuary.getDispersion(mobj);
            [Q,qIn,~] = Advection.getAdvectionFlow(mobj,'River');
            [Qtp,qtpIn,~] = Advection.getAdvectionFlow(mobj,'Qtp');
            %convert drifts to flow rates in m3/s
            qsIn = QSIn./cE/y2s;
            qsOut = QSOut./cE/y2s;
            dqIn = dExt+(kCeI.*qIn)+qtpIn+qsIn;
            qs = diag(1./cE)*QS/y2s;
            Qs = advMatrix(obj,qs,qsOut,nele);
            DQ = D+Q+Qtp+Qs;
            
            tol = cES/100; diff = 1; count = 0; conc0 = cES;
            while diff>tol && count<10  
                %iterate to find Qs based on convergence of concentration
                [conc,ok] = ASM_model.asmitaConcentrations(mobj,DQ,dqIn);
                if ok<1, qs = NaN; return; end %concentration matrix illconditioned
                qsOut = QSOut./conc/y2s;
                qs = diag(1./conc)*QS/y2s;
                Qs = advMatrix(obj,qs,qsOut,nele);
                DQ = D+Q+Qtp+Qs;
                diff = sum(abs(conc0-conc).*(qsOut~=0));
                conc0 = conc;
                count = count+1;
            end
        end
 
%%
        function ok = checkMassBalance(obj,AdvType)
            %check mass balance of advective flows
            obj = setAdvectionType(obj,AdvType);
            ok = checkAdvMassBalance(obj);
        end

%%  
        function initGraph = getInitialFlowGraph(obj,mobj,AdvType)
            %return a graph of the initial flow settings for AdvType
            [qin,qout,q] = getAdvectionProps(obj,AdvType);
            inoutxt = {'Downdrift','Updrift'};
            eleobj  = getClassObj(mobj,'Inputs','Element');
            nodetxt = setnodetext(eleobj,inoutxt);
            initGraph = matrix2graph(q,qin,qout,nodetxt);
        end

%%
        function [R,initGraph] = getFlowRatio(obj,mobj,AdvType)
            %find the ratio of the current flow rates to the initial values
            %R is struct of ratio (q/q0), diff (difference, q-q0) and 
            %diffratio (difference/initial value, dq/q0)                          
            drifts = obj.DriftGraph.Edges.Weight;
            initGraph = getInitialFlowGraph(obj,mobj,AdvType);
            initdrifts = initGraph.Edges.Weight;
            R.ratio = drifts./initdrifts;        
            R.diff = drifts-initdrifts;
            R.diffratio = abs(R.diff)./initdrifts;
        end

    end
%%
    methods (Static,Hidden)
        function setNewAdvection(mobj)
            %initialise an empty instance of Advection (used in asm_oo2mui)
            obj = Advection;
            setClassObj(mobj,'Inputs','Advection',obj);
        end        
    end
%%
%--------------------------------------------------------------------------
%       PRIVATE UTILITIES
%--------------------------------------------------------------------------      
    methods (Access=private)
        function [FlowIn,FlowOut,IntFlows] = getAdvectionProps(obj,AdvType)
            %return initial conditions of specified flow type as userdata
            switch AdvType
                case 'River'
                    FlowIn = obj.RiverIn;
                    FlowOut = obj.RiverOut;
                    IntFlows = obj.RiverFlows;
                case 'Drift'
                    FlowIn = obj.DriftIn;
                    FlowOut = obj.DriftOut;
                    IntFlows = obj.DriftFlows;
                case 'Qtp'
                    FlowIn = obj.TidalPumpingIn;
                    FlowOut = obj.TidalPumpingOut;
                    IntFlows = obj.TidalPumpingFlows;
            end
            %create a matrix comprising the internal flows and any inputs
            %and outputs, with dimensions [nele+2 x nele+2]
%             userdata = zeros(size(IntFlows));
%             userdata(end,2:end-1) = FlowIn;
%             userdata(2:end-1,1) = FlowOut;
%             userdata(2:end-1,2:end-1) = IntFlows;
        end

%%
        function obj = setAdvectionProps(obj,AdvType,idx)
            %assign the internal advection property values to the
            %specified set of flow type properties (River, Drift, Qtp)
            idin = 2; idout = 1;
            switch AdvType
                case 'River'
                    obj.RiverIn(idx,:) = [];
                    obj.RiverIn(:,idin) = obj.ExternalAdvIn;
                    obj.RiverOut(idx,:) = [];
                    obj.RiverOut(:,idout) = obj.ExternalAdvOut;
                    obj.RiverFlows = obj.InternalAdv;                    
                case 'Drift'
                    obj.DriftIn(idx,:) = [];
                    obj.DriftIn(:,idin) = obj.ExternalAdvIn;
                    obj.DriftOut(idx,:) = [];
                    obj.DriftOut(:,idout) = obj.ExternalAdvOut;
                    obj.DriftFlows = obj.InternalAdv;                    
                case 'Qtp'
                    idin = 1; idout = 2;
                    obj.TidalIn(idx,:) = [];
                    obj.TidalPumpingIn(:,idin) = obj.ExternalAdvIn;
                    obj.TidalOut(idx,:) = [];
                    obj.TidalPumpingOut(:,idout) = obj.ExternalAdvOut;
                    obj.TidalPumpingFlows = obj.InternalAdv;
            end
        end
%%
        function obj = addAdvectionProps(obj,AdvType)
            %assign the internal advection property values to the
            %specified set of flow type properties (River, Drift, Qtp)
            newlen = length(obj.ExternalAdvIn);
            AdvIn = zeros(newlen,2); AdvOut = AdvIn;
            if strcmp(AdvType,'Qtp')
                AdvIn(:,1) = obj.ExternalAdvIn;
                AdvOut(:,2) = obj.ExternalAdvOut;
            else
                AdvIn(:,2) = obj.ExternalAdvIn;
                AdvOut(:,1) = obj.ExternalAdvOut;
            end

            switch AdvType
                case 'River'
                    obj.RiverIn = AdvIn;
                    obj.RiverOut = AdvOut;
                    obj.RiverFlows = obj.InternalAdv;                    
                case 'Drift'
                    obj.DriftIn = AdvIn;
                    obj.DriftOut = AdvOut;
                    obj.DriftFlows = obj.InternalAdv;                    
                case 'Qtp'
                    obj.TidalPumpingIn = AdvIn;
                    obj.TidalPumpingOut = AdvOut;
                    obj.TidalPumpingFlows = obj.InternalAdv;
            end
        end
%%
        function addSource(obj,mobj,AdvType,idinput)
            %add a new advection source based on AdvType
            switch AdvType
                case 'River'
                    River.addRiver(mobj,idinput,...
                             obj.ExternalAdvIn(idinput,2));
                case 'Drift'
                    Drift.addDrift(mobj,idinput,...
                             obj.ExternalAdvIn(idinput,2));
            end
        end

%%                
        function obj = setAdvectionType(obj,AdvType)
            %assign the specified advection field to the internal
            %advection properties
            idin = 2; idout = 1;
            switch AdvType
                case 'River'
                    FlowIn = obj.RiverIn;
                    FlowOut = obj.RiverOut;
                    IntFlows = obj.RiverFlows;
                case 'Drift'
                    FlowIn = obj.DriftIn;
                    FlowOut = obj.DriftOut;
                    IntFlows = obj.DriftFlows;
                case 'Qtp'
                    FlowIn = obj.TidalPumpingIn;
                    FlowOut = obj.TidalPumpingOut;
                    IntFlows = obj.TidalPumpingFlows;
                    idin = 1; idout = 2;
            end
            %assign to generic input. Only need the source in and sink out
            obj.ExternalAdvIn = FlowIn(:,idin);   %column 2 is the source input
            obj.ExternalAdvOut = FlowOut(:,idout);%column 1 is the sink output
            obj.InternalAdv = IntFlows; 
        end
       
%%
        function clearAdvection(obj,AdvType)
            %clear advection properties based on AdvType
            switch AdvType
                case 'River'
                    obj.RiverIn = 0;
                    obj.RiverOut = 0;
                    obj.RiverFlows = [];
                    obj.RiverGraph = [];
                case 'Drift'
                    obj.DriftIn = 0;
                    obj.DriftOut = 0;
                    obj.DriftFlows = [];
                    obj.DrifGraph = [];
                case 'Qtp'
                    obj.TidalPumpingIn = 0;
                    obj.TidalPumpingOut = 0;
                    obj.TidalPumpingFlows = [];
                    obj.QtpGraph = [];
            end
             setClassObj(mobj,'Inputs','Advection',obj);
        end
        
%%
        function Q = advMatrix(~,q,qOut,nele)
            %set up advection flow matrix
            if isempty(q), Q = 0; return; end
            Q = -q;            
            for j=1:nele
                Q(j,j) = sum(q(j,:))+qOut(j);
            end
            % Use transpose of Q to match exchanges against concentrations
            Q = Q';
        end

%%
        function g_landward = tidalpumpingGraph(obj,mobj,qtp,qtp0)
            %extract the tidal pumping graph based on the tp discharges at the start
            %and end of each reach
            extype = mobj.GeoType(mobj.EXtypes);
            %invert the RiverGraph to get a landward version
            g_landward = inverse_graph(obj.RiverGraph);
            RchIDs = g_landward.Nodes.EleID;
            RchType = g_landward.Nodes.Type;
            RchName = g_landward.Nodes.Name;
            %to assign the qtp flows use the edge weights
            flowNodes = g_landward.Edges.EndNodes;
            flowWeight = g_landward.Edges.Weight;
            qtpWeight = zeros(size(flowWeight));

            isexternal = ismatch(RchType,extype);
            if any(isexternal)
                qtp = [qtp;0]; qtp0 = [0;qtp0];
            end
            %handle all the nodes linked to the outside
            outerNodes = successors(g_landward,1);            
            for i=1:length(outerNodes)  
                isextnode = strcmp(RchType(outerNodes(i)),extype);
                if any(isextnode)
                    %outerNode i is an external element type. Find what it
                    %links to in the network
                    chFlowID = successors(g_landward,outerNodes(i));                   
                else
                    %outerNode is a reach element linked to the outside
                    chFlowID = outerNodes(i);
                end
                chRchID = RchIDs(chFlowID);
                %need to remove external elements to get correct ID
%                 RchIDs(i) = [];
%                 innerRchs = RchIDs(2:end-1);
%                 chRchID = innerRchs==chRchID;
                qtpWeight(i) = qtp0(chRchID); %input to reaches from outside
            end
            
            innerNodes = find(qtpWeight==0);
            for i=1:length(innerNodes)  
                idx = innerNodes(i);      %iterate over list of innerNodes
                nodeName = flowNodes(innerNodes(i),1);
                idf = find(strcmp(RchName,nodeName));  %id of node i within full graph
                idr = find(strcmp(RchName(2:end-1),nodeName)); %Element Id (no offset)
                branchIDs = successors(g_landward,idf); %upstream branches connected to node i 
                if length(branchIDs)>1
                    %handle branches based on river flows
                    edgeidx = findedge(g_landward,idf,branchIDs);
                    totalFlow = sum(flowWeight(edgeidx)); %total in flow from connected branches
                    inFlow = flowWeight(idx); %in flow from ith branch 
                    qtpWeight(idx) = floor(qtp(idr)*inFlow/totalFlow);
                else
                    preBranchIDs = inedges(g_landward,idf);
                    if qtpWeight(preBranchIDs)<qtp(idr)
                        %trap branch with flow increaseing upstream - set 
                        %to downstream value
                        qtpWeight(idx) = sum(qtpWeight(preBranchIDs));
                    else
                        %direct connection
                        qtpWeight(idx) = qtp(idr);
                    end
                end                
            end
            
            g_landward.Edges.Weight = qtpWeight;
        end
%%
        function ok = checkAdvMassBalance(obj)
            %checks mass balance of properties loaded into internal
            %advection properties that depend on Advection Type
            %these need to be loaded before this function is called
            ok = 1;
            q = obj.InternalAdv;
            qIn = obj.ExternalAdvIn;
            qOut = obj.ExternalAdvOut;
            nele = length(qIn);
            for i=1:nele
                sumq = -sum(q(i,:))+ sum(q(:,i))+ qIn(i) - qOut(i);
                if sumq~=0
                    ok = 0;
                end
            end
        end   

    end
end
