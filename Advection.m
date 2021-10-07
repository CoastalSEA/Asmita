classdef Advection < handle
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
        RiverIn = 0            %vector of external flows into estuary
        RiverOut = 0           %vector of external flows out of estuary
        DriftFlows = []        %matrix of littoral drift flows, qs (m3/s)
        DriftIn = 0            %vector of littoral flows into system
        DriftOut = 0           %vector of littoral flows out of system
    end
    
    properties (Transient)  %properties that are time varying
        TidalPumpingFlows = [] %tidal pumping matrix, qtp (m3/s)
        TidalPumpingIn = 0     %vector of tidal pumping flows into estuary
        TidalPumpingOut = 0    %vector of tidal pumping flows out of estuary
        RiverGraph             %handle to graph of river flows (advection)
        DriftGraph             %handle to graph of drift flows (advection)
        QtpGraph               %handle to graph of tidal pumping flows (advection)  
    end
    
    properties (Access=private, Transient)
        InternalAdv = []         %advection matrix, q (m3/s)
        ExternalAdvIn = 0        %matrix of external flows into estuary
        ExternalAdvOut = 0       %matrix of external flows out of estuary        
    end
    
    methods (Access=private)
        function obj = Advection
        end
    end
    
    methods (Static)
%--------------------------------------------------------------------------
%       SETUP
%--------------------------------------------------------------------------
        function obj = getAdvection
            %test whether Model exists. Only allow single instance
            %the following code forces a singleton class
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = Advection;
                %constructor code for advection class
            end
            obj = localObj;
        end
                      
%%        
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
            elseif isempty(rivobj) && isempty(dftobj)
                warndlg('No river or drift inputs have beend defined')
                return; 
            end
            nele = length(eleobj);
            nelp = nele+2;
         
            %create or get data and build matrix of exchanges
            userdata = zeros(nelp,nelp);
            obj = Advection.setAdvectionType(mobj,AdvType);
            if ~isempty(obj.InternalAdv)                
                userdata(end,2:end-1) = obj.ExternalAdvIn;
                userdata(2:end-1,1) = obj.ExternalAdvOut;
                userdata(2:end-1,2:end-1) = obj.InternalAdv;
            else
                obj.ExternalAdvIn = zeros(nele,1);
            end
            %
            InputEle = []; flow = 0;
            switch AdvType
                case 'River'
                    if ~isempty(rivobj)
                        h_flow = rivobj;
                        flowVar = 'RiverFlow';
                        [InputEle,okEle,msg,flow] = River.getSourceProps(mobj);     
                    end
                    qunits = '(m3/s)';
                    inoutxt = {'0: Outside';'+: Rivers'};
                case 'Drift'
                    if ~isempty(dftobj)
                        h_flow = dftobj;
                        flowVar = 'DriftRate';
                        [InputEle,okEle,msg,flow] = Drift.getSourceProps(mobj); 
                    end
                    qunits = '(m3/year)';
                    inoutxt = {'>: Sink';'>: Source'};
            end
            MatrixInpEle = InputEle+1;         %offset due to 'outside'
            userdata(nelp,MatrixInpEle) = flow;

            prop = 'Advection';
            prompt = sprintf('Enter advection between elements %s\n(from row element to column element in direction of flow)',qunits);            
%             h_tab = Estuary.setMatrix(mobj,prop,prompt,inoutxt,userdata);
            matrixtable = setMatrix(eleobj,prop,prompt,inoutxt,userdata);
%             if isempty(matrixtable), return; end %user cancelled
            %write user defined date to variables
            preExternalAdvIn = obj.ExternalAdvIn;
            extadv = matrixtable{end,2:end-1};   %from outside
            obj.ExternalAdvIn = extadv';
            extadv = matrixtable{2:end-1,1};     %to outside
            obj.ExternalAdvOut = extadv;
            obj.InternalAdv = matrixtable{2:end-1,2:end-1};

            %check whether inputs have been changed and update
            idx = diff([preExternalAdvIn,obj.ExternalAdvIn],1,2);
            idinputs = find(idx~=0);
            if ~isempty(InputEle)     %if some source already exist
                idnewinput = setdiff(idinputs,InputEle);
                if ~isempty(idnewinput)
                    for i=1:length(idnewinput)   %add new sources                        
                        %check that source is assigned to an element that exists   
                        ok = okEle(idnewinput(i));
                        if ok==0 %catch changes to elements that are not valid
                            if obj.ExternalAdvIn(idnewinput(i))>0 %only report when a non-zero value has been input (ie not corrections back to zero)
                                msgtxt = sprintf(msg{1},idnewinput(i));
                                warndlg(msgtxt);
                            end
                        else
                            %add a new source
                            Advection.addSource(mobj,AdvType,idnewinput(i));
                            msgtxt = sprintf(msg{2},idnewinput(i));
                            warndlg(msgtxt);
                        end
                    end
                end
                for i=1:length(InputEle)
                    h_flow(i).(flowVar) =  obj.ExternalAdvIn(InputEle(i));
                end
                
            else                        %if no input source exist
                 if ~isempty(idinputs)
                     warndlg('Inputs have been added and need to be defined in Setup');
                     for i=1:length(idinputs)
                         Advection.addSource(mobj,AdvType,idinputs(i));
                     end
                 end
            end
            %assign updated flow field to the specified flow type
            obj = setAdvectionProps(obj,AdvType);
            obj = setAdvectionGraph(obj,mobj,AdvType);    
            setClassObj(mobj,'Inputs','Advection',obj);
        end

%%
        function addSource(mobj,AdvType,idinput)
            %add a new advection source based on AdvType
            obj  = getClassObj(mobj,'Inputs','Advection');
            switch AdvType
                case 'River'
                    River.addRiver(mobj,idinput,...
                             obj.ExternalAdvIn(idinput));
                case 'Drift'
                    Drift.addDrift(mobj,idinput,...
                             obj.ExternalAdvIn(idinput));
            end
        end
        
%%       
        function addEleAdvection(mobj)
            %amend advection matrix to add an extra element
            obj  = getClassObj(mobj,'Inputs','Advection');
            if ~isempty(obj.RiverFlows)
                obj = Advection.setAdvectionType(mobj,'River');
                obj = updateAdvectionArrays(obj);
                obj = setAdvectionProps(obj,'River');
            end
            %
            if ~isempty(obj.DriftFlows)
                obj = Advection.setAdvectionType(mobj,'Drift');
                obj = updateAdvectionArrays(obj);
                obj = setAdvectionProps(obj,'Drift');
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
        function delEleAdvection(mobj,idx)
            %delete advection properties for element(idx)
            obj  = getClassObj(mobj,'Inputs','Advection');
            if ~isempty(obj.RiverFlows)
                obj = Advection.setAdvectionType(mobj,'River');
                obj = deleteAnAdvection(obj,idx);
                obj = setAdvectionProps(obj,'River');
            end
            %
            if ~isempty(obj.DriftFlows)
                obj = Advection.setAdvectionType(mobj,'Drift');
                obj = deleteAnAdvection(obj,idx);
                obj = setAdvectionProps(obj,'Drift');
            end
            %
            %assign updated flow field to the specified flow type
            setClassObj(mobj,'Inputs','Advection',obj);
            %--------------------------------------------------------------
            function obj = deleteAnAdvection(obj,idx)
                obj.ExternalAdvIn(idx) = [];
                obj.ExternalAdvOut(idx) = [];
                obj.InternalAdv(idx,:) = [];
                obj.InternalAdv(:,idx) = [];
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
                nodetxt = Estuary.setGraphVariables(mobj,'Outside','River(s)');
                flowGraph = Advection.setFlowGraph(mobj,Adv,AdvIn,AdvOut,...
                                                             nodetxt);
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
                    initialiseFlow(rivobj);
                end                                         
                nlabel = strcat(num2str(flowGraph.Nodes.EleID),...
                    '-',flowGraph.Nodes.Name);
                %check mass balance of advective flows
                ok = Advection.checkMassBalance(mobj,'River');
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
                nodetxt = Estuary.setGraphVariables(mobj,'Downdrift','Updrift');
                driftGraph = Advection.setFlowGraph(mobj,Qs,QsIn,QsOut,...
                                                            nodetxt);
                dftobj = getClassObj(mobj,'Inputs','Drift');                                    
                j = zeros(length(dftobj),1); 
                for i=1:length(dftobj) %find souces defined by a timeseries
                   j(i) = ~isempty(dftobj(i).DriftTSC); 
                end    
                %
                if sum(j)>0  %update if timeseries input
                    obj.DriftGraph = driftGraph; %update the model instance
                    setClassObj(mobj,'Inputs','Advection',obj);
                    %
                    rnpobj = getClassObj(mobj,'Inputs','RunProperties');
                    Advection.updateAdvectionGraphs(mobj,rnpobj.StartYear);
                    obj  = getClassObj(mobj,'Inputs','Advection');
                    driftGraph = obj.DriftGraph;
                else
                    initialiseFlow(dftobj);
                end                                                        %                                       
                nlabel = strcat(num2str(driftGraph.Nodes.EleID),...
                    '-',driftGraph.Nodes.Name);
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
                %intitialise transient properties
                AsmitaModel.intialiseModelParameters(mobj);
                %Advection.setTidalPumping(mobj);     %initialise qtp                    
                if isempty(obj.TidalPumpingFlows) %if still not defined
                    msg = sprintf('Check parameters for tidal pumping have been defined');
                    nlabel = {msg};
                    return;
                end
                Qtp = obj.TidalPumpingFlows;
                QtpIn = obj.TidalPumpingIn;
                QtpOut = obj.TidalPumpingOut;
                nodetxt = Estuary.setGraphVariables(mobj,'River(s)','Outside');
                qtpGraph = Advection.setFlowGraph(mobj,Qtp,QtpIn,QtpOut,...
                                                                nodetxt);
                nlabel = strcat(num2str(qtpGraph.Nodes.EleID),...
                    '-',qtpGraph.Nodes.Name);
                %check mass balance of tidal pumping flows
%                 ok = Advection.checkMassBalance(mobj,'Qtp');
%                 if ok<1
%                     warndlg('Mass balance fails for defined Tidal Pumping flows');
%                 end
            end
        end
       
%%
%--------------------------------------------------------------------------
%       RUNTIME
%--------------------------------------------------------------------------
        function [Q,qIn,qOut] = getAdvectionFlow(mobj,AdvType)
            %get the advection properties as an element matrix
            AdvGraphName = strcat(AdvType,'Graph');     
            obj  = getClassObj(mobj,'Inputs','Advection');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if isempty(obj) || isempty(obj.(AdvGraphName)) || isempty(obj.(AdvGraphName).Nodes)
                Q = 0; qIn = 0; qOut = 0;
                return;
            end
            %get flow properties using graph if dynamic or properties if not
            if rncobj.IncDynamicElements
                AdvGraph = obj.(AdvGraphName);
                varmatrix = Advection.getSubVarMatrix(mobj,AdvGraph);
                qIn = varmatrix(end,2:end-1)';
                qOut = varmatrix(2:end-1,1);
                q = varmatrix(2:end-1,2:end-1);
            else 
                obj = Advection.setAdvectionType(mobj,AdvType);
                qIn = obj.ExternalAdvIn;
                qOut = obj.ExternalAdvOut;
                q = obj.InternalAdv;
            end
            
            if strcmp(AdvType,'Drift')
                [q,qIn,qOut] = Advection.getDriftFlow(mobj,q,qIn,qOut);
            end
            
            %set up flow matrix
            Q = -q;
            nele = length(getClassObj(mobj,'Inputs','Element'));
            for j=1:nele
                Q(j,j) = sum(q(j,:))+qOut(j);
            end
            % Use transpose of Q to match exchanges against concentrations
            Q = Q'; 
        end 
        
%%
        function [Qs,qsIn,qsOut] = getDriftFlow(mobj,QS,QSIn,QSOut)            
            %adjust littoral drift rates to equivalent flow rates
            eleobj = getClassObj(mobj,'Inputs','Element');
            estobj = getClassObj(mobj,'Inputs','Estuary');
            cE = getEleProp(eleobj,'EqConcentration');
            kCeI = River.getRiverProp(mobj,'tsRiverConc')./cE; 
            y2s = mobj.Constants.y2s;
            cES = estobj.EqConcCoarse;
            [D,dExt] = Estuary.getDispersion(mobj);
            [Q,qIn,~] = Advection.getAdvectionFlow(mobj,'River');
            [Qtp,qtpIn,~] = Advection.getAdvectionFlow(mobj,'Qtp');
            %convert drifts to flow rates in m3/s
            qsIn = QSIn./cE/y2s;
            dqIn = dExt+kCeI.*qIn+qtpIn+qsIn;
            Qs = diag(1./cE)*QS/y2s;
            DQ = D+Q+Qtp+Qs;
            tol = cES/100; diff = 1; count = 0; conc0 = cES;
            while diff>tol && count<10  
                %iterate to find Qs based on convergence of concentration
                conc = ASM_model.asmitaConcentrations(mobj,DQ,dqIn);
                qsOut = QSOut./conc/y2s;
                Qs = diag(1./conc)*QS/y2s;
                DQ = D+Q+Qtp+Qs;
                diff = sum(abs(conc0-conc).*(qsOut~=0));
                conc0 = conc;
                count = count+1;
            end
        end
        
%% 
        function updateAdvectionGraphs(mobj,tsyear)
            %initialise and update dynamic graphs based on modified conditions (if any)            
            obj  = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj)
                return;
            end
            rivobj  = getClassObj(mobj,'Inputs','River');
            %only check need to update if river is included
            if ~isempty(obj.RiverGraph) && ~isempty(rivobj)
                obj.RiverGraph = getflowpath(rivobj.River);                
            end
            dftobj  = getClassObj(mobj,'Inputs','Drift');
            %only check need to update if drift is included
            if ~isempty(obj.DriftGraph) && ~isempty(dftobj)
                obj.DriftGraph = getflowpath(dftobj.Drift);
            end

             setClassObj(mobj,'Inputs','Advection',obj);
            %--------------------------------------------------------------
            function h_flowpath = getflowpath(gobj)
                %internal function to get relevant flowgraph
                if isa(gobj,'River')
                    GraphName = 'RiverGraph';
                    TSname = 'RiverTSC';
                elseif isa(gobj,'Drift')
                    GraphName = 'DriftGraph';
                    TSname = 'DriftTSC';
                else
                    return;
                end
                h_flowpath = obj.(GraphName); %existing definition
                j = zeros(length(gobj),1); 
                for is=1:length(gobj) %find sources defined by a timeseries
                   j(is) = ~isempty(gobj(is).(TSname)); 
                end                
                if sum(j)>0  %timeseries present so update flowgraph
                    for is=1:length(j) %loop to consider each source                        
                        if j(is)>0
                            %interpolate time series to get new input value
                            if isa(gobj,'Drift') 
                                [inFlow,inEleID] = Drift.getDriftTSprop(mobj,...
                                                    is,tsyear);
                            else
                                [inFlow,inEleID] = River.getRiverTSprop(mobj,...
                                                    is,tsyear);
                            end
                            %
                            if length(inFlow)>1 %drifts vary along path
                                h_flowpath = newDriftPath(h_flowpath,...
                                                        inFlow,inEleID);
                            else %source value dictates flow along path
                                h_flowpath = newFlowPath(h_flowpath,...
                                                        inFlow,inEleID);  
                            end
                        end
                    end
                end 
            end
            %--------------------------------------------------------------
            function h_flowpath = newFlowPath(h_flowpath,newInFlow,inEleID)
                %internal function to reassign new flows to flowgraph
                %used for all river flows and drifts when source value 
                %dictates flow along path
                flowpathID = h_flowpath.Nodes.EleID;
                %find flow input
                inid = find(flowpathID==inEleID);
                newEdgeWeight = h_flowpath.Edges.Weight;
                %define source input as first link
                %naming convention is *ID for elements and *id for h_flowpath index
                %
                %id that flow/drift is coming from (source) 
                upstreamids = predecessors(h_flowpath,inid);
                outid = upstreamids((flowpathID(upstreamids)==0));
                ide = findedge(h_flowpath,outid,inid);
                newEdgeWeight(ide) = newInFlow;
                while ~isempty(inid)
                    nextin = [];
                    cumOutFlow = zeros(length(inid),1);
                    for j = 1:length(inid) %handle branching
                        %handle flowpath input
                        ide = findedge(h_flowpath,outid,inid(j));
                        oldInFlow = h_flowpath.Edges.Weight(ide);                        
                        %find other inputs to the same node.
                        %naming convention *Input = sum(*InFlow)
                        upstreamids = predecessors(h_flowpath,inid(j));
                        idf = findedge(h_flowpath,upstreamids,inid(j));
                        oldInput = sum(h_flowpath.Edges.Weight(idf)); %total input
                        if length(ide)>1 %handle two new inputs (diverge/recvonverge case)
                            newInput = oldInput-sum(oldInFlow)+sum(newInFlow);
                        else
                            newInput = oldInput-oldInFlow+newInFlow(j);
                        end
                        %find outflowing edges and update pro-rata based on
                        %previous flow distribution. Naming convention is
                        %*Output = sum(*OutFlow)
                        downstreamids = successors(h_flowpath,inid(j));
                        if ~isempty(downstreamids)
                            ido = findedge(h_flowpath,inid(j),downstreamids);
                            oldOutFlow = h_flowpath.Edges.Weight(ido);
                            oldOutput = sum(oldOutFlow);
                            newOutFlow = zeros(length(downstreamids),1);
                            %weight allows drift to vary based on initial 
                            %definition in advection matrix. Does not
                            %affect river flows because mass balance imposed
                            weight = oldOutput/oldInput; 
                            for i=1:length(downstreamids)
                                newOutFlow(i) = oldOutFlow(i)/oldOutput*newInput*weight;
                                newEdgeWeight(ido(i)) = newOutFlow(i);
                            end
                            nextin = [nextin,downstreamids]; %#ok<AGROW>
                            cumOutFlow(j) = sum(newOutFlow);
                        end
                    end
                    outid = inid(flowpathID(inid)~=0);
                    inid = unique(nextin); 
                    %handle special case of links that diverge and reconverge
                    if length(outid)>1 
                        newInFlow = cumOutFlow;
                    else
                        newInFlow = newOutFlow;
                    end                                       
                end
                h_flowpath.Edges.Weight = newEdgeWeight;
            end
            %--------------------------------------------------------------
            function h_flowpath = newDriftPath(h_flowpath,newDrift,inEleID)
                %internal function to reassign new flows to flowgraph
                %when there are multiple sources. No checks on mass balance
                flowpathID = h_flowpath.Nodes.EleID;
                newEdgeWeight = h_flowpath.Edges.Weight;
                %find Edge for each dynamic drift input ans assign update
                for j=1:length(newDrift)
                    if inEleID(j,1)==0                        
                        outID = find(flowpathID==inEleID(j,2));
                        upstreamids = predecessors(h_flowpath,outID);
                        inID = upstreamids((flowpathID(upstreamids)==0));
                    elseif inEleID(j,2)==0
                        inID = find(flowpathID==inEleID(j,1));
                        dnstreamids = successors(h_flowpath,inID);
                        outID = dnstreamids((flowpathID(dnstreamids)==0));
                    else
                        inID = find(flowpathID==inEleID(j,1));
                        outID = find(flowpathID==inEleID(j,2));
                    end
                    ide = findedge(h_flowpath,inID,outID);
                    newEdgeWeight(ide) = newDrift(j);
                end
                h_flowpath.Edges.Weight = newEdgeWeight;
            end
            %
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
                if isempty(obj.RiverGraph)
                    %intitialise RiverGraph
                    obj.RiverGraph = Advection.initialiseRiverGraph(mobj);
                end                
            end

            %allocate fluxes to the correct postion in matrix
            %to find correct node use a landward version of FlowGraph
            h_landward = Reach.LandwardFlow(obj.RiverGraph);
            FlowPathIDs = h_landward.Nodes.EleID;
            EleType = h_landward.Nodes.Type;
            EleName = h_landward.Nodes.Name;
            %delFlowID = find(strcmp(EleType,'Delta'));
            %to access reach properties use the landward channel paths
            h_channels = Reach.LandwardPath(mobj);  %uses channels in current version of ReachGraph            
            reachChannelIDs = h_channels.Nodes.EleID;
            reachChannelName = h_channels.Nodes.Name;
            
            %get reach based tidal pumping discharge
            [qtp,qtp0] = Advection.getTidalPumpingDischarge(mobj,reachChannelIDs);   
            
            %to assign the qtp flows use the edge nodes
            flowNodes = h_landward.Edges.EndNodes;
            flowWeight = h_landward.Edges.Weight;
            qtpWeight = zeros(size(flowWeight));
            qtpOuter = qtpWeight;
            %handle all the nodes linked to the outside
            outerNodes = successors(h_landward,1);
            edgeNodes = strcmp(flowNodes(:,1),'Outside');
            outerEles = flowNodes(edgeNodes,2);
            for i=1:length(outerNodes)
                if strcmp(EleType(outerNodes(i)),'Delta')
                    %delEleID  = FlowPathIDs(outerlinks(i)); %used as check
                    chFlowID = successors(h_landward,outerNodes(i));
                else
                    chFlowID = outerNodes(i);
                end
                chEleID = FlowPathIDs(chFlowID);
                chRchID = reachChannelIDs==chEleID;  
                qtpWeight(i) = qtp0(chRchID);
                qtpOuter(i) = 1;
                if strcmp(EleType(outerNodes(i)),'Delta')
                    idx = strcmp(flowNodes(:,1),outerEles(i));                    
                    qtpWeight(idx) = qtp0(chRchID);
                    qtpOuter(idx) = 1;
                end
            end

            innerNodes = find(qtpOuter<1);
            for i=1:length(innerNodes)
                idx = innerNodes(i);
                nodeName = flowNodes(innerNodes(i),1);
                idf = find(strcmp(EleName,nodeName));
                idr = find(strcmp(reachChannelName,nodeName));
                branchIDs = successors(h_landward,idf);
                if length(branchIDs)>1
                    %handle branches based on river flows
                    edgeidx = findedge(h_landward,idf,branchIDs);
                    totalFlow = sum(flowWeight(edgeidx));
                    inFlow = flowWeight(idx);
                    %check whether this should be (1-ratio)********************
                    qtpWeight(idx) = floor(qtp(idr)*inFlow/totalFlow);
                else
                    qtpWeight(idx) = qtp(idr);
                end
            end
            
            h_landward.Edges.Weight = qtpWeight;
            varmatrix = Advection.getSubVarMatrix(mobj,h_landward);
            qtpIn = varmatrix(end,2:end-1)';
            qtpOut = varmatrix(2:end-1,1);
            qtpFlows = varmatrix(2:end-1,2:end-1);
            %assign instance values to internal properties
            obj.InternalAdv = qtpFlows;    
            obj.ExternalAdvIn = qtpIn;      
            obj.ExternalAdvOut = qtpOut;
            %assign updated flow field to the tidal pumping properties
            obj = setAdvectionProps(obj,'Qtp');
            obj = setAdvectionGraph(obj,mobj,'Qtp'); 
            
            setClassObj(mobj,'Inputs','Advection',obj);
        end 
        
%%
%--------------------------------------------------------------------------
%       UTILITIES
%--------------------------------------------------------------------------
        function g = setFlowGraph(mobj,adv,advIn,advOut,nodetxt)
            %assign the flow graph to handle so that it can be used elsewhere
            %uses flow matrix, adv, with inputs, advin, and outputs,
            %advout, to define flow paths. nodetxt provided by calling function              
            eleobj  = getClassObj(mobj,'Inputs','Element');
            nele = length(eleobj);
            userdata(nele+2,:) = zeros(1,nele+1);
            userdata(:,nele+2) = zeros(nele+2,1);
            userdata(end,2:end-1) = advIn;
            userdata(2:end-1,1) = advOut;
            userdata(2:end-1,2:end-1) = adv;

            if sum(sum(userdata))==0
                userdata = 0;
%                 nname = {'No Flows Defined'};
            end
            g = digraph(userdata);
            %now remove elements that are not connected
            [in,out] = findedge(g);
            idx = unique([in,out]);
            g = subgraph(g,idx);
            g.Nodes.EleID = nodetxt.nid(idx);
            g.Nodes.Type = nodetxt.ntype(idx);
            %node names must be unique
            nname = Estuary.checkUniqueNames(nodetxt.nname);
            g.Nodes.Name = nname(idx);
        end

%%
        function varmatrix = getSubVarMatrix(mobj,vargraph)  
            %extract components of graph to define advection   
            eleobj  = getClassObj(mobj,'Inputs','Element');
            nele = length(eleobj);
            [s,t] = findedge(vargraph);
            sEleID = vargraph.Nodes.EleID(s);
            tEleID = vargraph.Nodes.EleID(t);
            EleID = getEleProp(eleobj,'EleID');
            sEleIdx = ones(size(sEleID)); tEleIdx = sEleIdx;
            for i=1:length(sEleID)
                if sEleID(i)>0
                    sEleIdx(i) = find(EleID==sEleID(i))+1;
                else
                    sEleIdx(i) = nele+2;
                end
                if tEleID(i)>0
                    tEleIdx(i) = find(EleID==tEleID(i))+1;
                end
            end
            varmatrix = sparse(sEleIdx,tEleIdx,vargraph.Edges.Weight,nele+2,nele+2);
            varmatrix = full(varmatrix);
        end
        
%%                
        function obj = setAdvectionType(mobj,AdvType)
            %define or update advection properties
            obj  = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj)
                obj = Advection.getAdvection;
            end
            %assign the specified advection field to the internal
            %advection properties
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
            obj.ExternalAdvIn = FlowIn;
            obj.ExternalAdvOut = FlowOut;
            obj.InternalAdv = IntFlows; 
        end
  
%%
        function [qtp,qtp0] = getTidalPumpingDischarge(mobj,rchChIDs)
            %compute tidal pumping for each reach with a river flow based on tide and river properties
            wlvobj  = getClassObj(mobj,'Inputs','WaterLevels');
            aM4  = wlvobj.M4amplitude; %amplitude of M4 tidal constituent (m)        
            phM4 = wlvobj.M4phase;     %phase of M4 tidal constituent (m)            
            zHW0 = wlvobj.HWaterLevel; %High Water Level at mouth at time t         
            zLW0 = wlvobj.LWaterLevel; %Low Water Level at mouth at time t
            Tp   = wlvobj.TidalPeriod; %tidal period (hrs)            

            eleobj = getClassObj(mobj,'Inputs','Advection');
            eletype = getEleProp(eleobj,'transEleType');
            ich   = strcmp(eletype,'Channel');
            eLe = getEleProp(eleobj,'Length');
            ws = getEleProp(eleobj,'transVertExch');
            cero = getEleProp(eleobj,'SedMobility');
            
            estobj = getClassObj(mobj,'Inputs',Estuary');
            eLw  = estobj.WidthELength;%estuary width e-folding length (m)
            
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
            avCSA = Reach.getReachProp(mobj,'UpstreamCSA');
            %derive river flow speed from river discharge and local(upstream)csa
            uriv = Reach.getReachProp(mobj,'RiverFlow');
            %
            % tpump = (3/2.cero.U.(aeff.k.x/H-(a2/a).sin(phM4)+1/4.Tc.U/Lw)-uriv*(3/2.cero-1)).Axsec  [eq.1]
            % U  = sig*Lw*a/H; Tc = H/ws;
            % cero is defined as alpha' in the source paper and manual
            % compute tidal pumping flux per unit area. Re-write eqn.1 per unit area as:
            % tpump'=1.5*cero*sig*eLw*amp/Hav*(aef*kno*xi/Hav-M4M*sin(phM4)+0.25*sig*amp/ws)-uriv*(3/2.cero-1)
            term1 = 3/2*cero(ich)*U;
            term2 = aef.*kno.*xi/Hav;
            term3 = M4M*sin(phM4*pi/180); 
            term4 = 0.25*sigma*amp./ws(ich);
            term5 = uriv.*(3/2*cero(ich)-1);
            qtp = floor((term1.*(term2-term3+term4)-term5).*avCSA);
            %area and tidal pumping at the mouth of each inlet            
            Am0 = avCSA.*exp(eLe(rchChIDs)/eLw);
            qtp0 = floor((term1.*(0-term3+term4)-term5).*Am0); 
        end

%%
        function clearAdvection(mobj,AdvType)
            %clear advection properties based on AdvType
            obj  = getClassObj(mobj,'Inputs','Advection');
            if isempty(obj), return; end
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
        function ok = checkMassBalance(mobj,AdvType)
            %check mass balance of advective flows
            obj = Advection.setAdvectionType(mobj,AdvType);
            ok = checkAdvMassBalance(obj);
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

%--------------------------------------------------------------------------
%       PRIVATE UTILITIES
%--------------------------------------------------------------------------      
    methods (Access=private)
        function userdata = getAdvectionProps(obj,AdvType)
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
            userdata = zeros(size(IntFlows));
            userdata(end,2:end-1) = FlowIn;
            userdata(2:end-1,1) = FlowOut;
            userdata(2:end-1,2:end-1) = IntFlows;
        end

%%
        function obj = setAdvectionProps(obj,AdvType)
            %assign the internal advection property values to the
            %specified set of flow type properties (River, Drift, Qtp)
            switch AdvType
                case 'River'
                    obj.RiverIn = obj.ExternalAdvIn;
                    obj.RiverOut = obj.ExternalAdvOut;
                    obj.RiverFlows = obj.InternalAdv;                    
                case 'Drift'
                    obj.DriftIn = obj.ExternalAdvIn;
                    obj.DriftOut = obj.ExternalAdvOut;
                    obj.DriftFlows = obj.InternalAdv;                    
                case 'Qtp'
                    obj.TidalPumpingIn = obj.ExternalAdvIn;
                    obj.TidalPumpingOut = obj.ExternalAdvOut;
                    obj.TidalPumpingFlows = obj.InternalAdv;
            end
        end
        
%%
        function obj = setAdvectionGraph(obj,mobj,AdvType)
            %initialise a flow graph based on AdvType and the current 
            %internal advection properties
            AdvGraph = strcat(AdvType,'Graph');
            switch AdvGraph
                case 'RiverGraph'
                    startname = 'Outside';
                    endname = 'River(s)';
                case 'DriftGraph'
                    startname = 'Updrift';
                    endname = 'Downdrift';
                case 'QtpGraph'
                    startname = 'River(s)';
                    endname = 'Outside';
            end
            Adv = obj.InternalAdv;
            AdvIn = obj.ExternalAdvIn;
            AdvOut = obj.ExternalAdvOut;
            nodetxt = Estuary.setGraphVariables(mobj,startname,endname);
            flowGraph = Advection.setFlowGraph(mobj,Adv,AdvIn,AdvOut,...
                                                                nodetxt);
            obj.(AdvGraph) = flowGraph;
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
