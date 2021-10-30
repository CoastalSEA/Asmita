classdef Reach < handle
    %Reaches are used to hold various gross properties that are used in the
    %model and results output. A reach is defined as a single channel 
    %element and all other types of element that are linked to it. 
    %
    %----------------------------------------------------------------------
    % AUTHOR
    % Ian Townend
    %
    % COPYRIGHT
    % CoastalSEA, (c) 2016
    %----------------------------------------------------------------------
    %  
    properties (Access=private, Transient)
        DispersionGraphID        %Elements within a reach as defined in DispersionGraph
    end
    
    properties (Transient)
        ReachChannelID      %Channel elements that define each reach
        ReachEleID          %Elements within a reach using Element ID
        HWvolume            %total volume of reach to high water
        LWvolume            %total volume to low water (usually just channel)
        HWarea              %total plan area of reach at high water
        LWarea              %total plan area of reach at low water
        HWlevel             %elevation of HW in reach
        LWlevel             %elevation of LW in reach
        MWlevel             %elevation of MWL in reach
        ReachPrism          %volume between high and low water in reach
        UpstreamPrism       %volume between high and low water of current 
                            %reach and all reaches upstream of current reach
        ReachLength         %length of reach (m)
        CumulativeLength    %distance to upstream end of element from mouth (excl.delta)
        ReachCSA            %cross-sectional area of reach (m2)
        UpstreamCSA         %cross-sectional area of upstream end of reach at mtl
        RiverFlow           %Velocity at upstream end of reach based on river discharge and CSA
    end
    
    methods (Access=private)
        function obj = Reach
        end
    end
    
%%    
    methods (Static)
        function setReach(mobj)
            %initialise reach and set up DispersionGraph and RiverGraph based on
            %initial conditions
            classname = 'Reach';
            obj = getClassObj(mobj,'Inputs',classname);            
            if ~isa(obj,'Reach')
                obj = Reach;
            end
            setClassObj(mobj,'Inputs','Reach',obj);
            Reach.setReachProps(mobj); %also called directly in AsmitaModel.InitTimeStep
            
            %assign initial reach id to elements for use in plotting   
            eleobj = getClassObj(mobj,'Inputs','Element');
            for i=1:length(obj)
                eleobj(obj(i).ReachChannelID).ReachID = i;
                for j=1:length(obj(i).ReachEleID)
                    eleobj(obj(i).ReachEleID(j)).ReachID = i;
                end
            end
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
%%
        function setReachProps(mobj)
            %assign volume and area properties to each reach
            %assumes that reach links can change so reassigns each time
            %the function is called using values in DispersionGraph and RiverGraph 
            obj = getClassObj(mobj,'Inputs','Reach');
            if isempty(obj)
                obj = Reach.setReach(mobj);
            end
            %additonal class instances needed
            eleobj = getClassObj(mobj,'Inputs','Element');
            estobj = getClassObj(mobj,'Inputs','Estuary');
            advobj = getClassObj(mobj,'Inputs','Advection');
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            cstobj = getClassObj(mobj,'Inputs','CSThydraulics');

            %element and water level properties
            V = getEleProp(eleobj,'MovingVolume');
            S = getEleProp(eleobj,'MovingSurfaceArea');
            n = getEleProp(eleobj,'TransportCoeff');
            eLe = getEleProp(eleobj,'Length');
            eleID = getEleProp(eleobj,'EleID');
            damp = getEleProp(eleobj,'TidalDamping');
            if isempty(V)
                warndlg('No Element Volumes defined');
                return
            end
            
            if isempty(advobj.RiverGraph)
                %advection network has not been initialised
                advobj.RiverGraph = Advection.initialiseRiverGraph(mobj);
            end
            
            %connectivity graphs and type + node ids
            h_network = estobj.DispersionGraph; %full network in landward direction
            eletype = h_network.Nodes.Type;
            reachEleID = h_network.Nodes.EleID;
            h_landward = type_sub_graph(h_network,'Channel'); %reaches in landward direction
%             h_landward = Reach.LandwardPath(mobj);  
            reachChannelID = h_landward.Nodes.EleID;

            for i=1:length(reachChannelID) %uses graph ids and element ids
                %set channel based properties of reach
                idx = eleID==reachChannelID(i);
                obj(i,1).LWvolume = V(idx);
                obj(i,1).LWarea = S(idx);                
                obj(i,1).ReachLength = eLe(idx);
            end
            
            %set reach length of each element and cumulative length from mouth
            setReachLength(obj,reachChannelID,h_landward)
            
            %set channel based properties of each reach
            setReachChannelProperties(obj,reachChannelID,reachEleID,...    
                                                         eletype,h_network);

            %set high and low water levels using output of CSTmodel if used
            [HWL,MWL,LWL,TR] = setHighLowWater(obj,reachChannelID,...
                                        rncobj,cstobj,advobj,wlvobj,damp);

            %set collective properties for each reach
            setCollectiveReachProperties(obj,reachChannelID,eleID,eletype,...
                                                    V,S,n,HWL,MWL,LWL,TR);

            %set path dependent upstream CSA (uses graph ids)
            setUpstreamCSA(obj,reachChannelID,h_landward);

            %river flow velocity at upstream end of element
            setRiverFlow(obj,reachChannelID,advobj);

            %cumulative upstream tidal prism (including reach prism of element)
            setCumTidalPrism(obj,reachChannelID,h_landward)

            setClassObj(mobj,'Inputs','Reach',obj);
        end

%% ------------------------------------------------------------------------
% functions called to access properties from external methods
%--------------------------------------------------------------------------      
        function prop = getReachProp(mobj,varname)
            %property varname returned as a reach array
            obj = getClassObj(mobj,'Inputs','Reach');
            if isempty(obj) || isempty(obj(1).(varname))
                %handle case of a channel only reach when ReachEleID is
                %empty but still need to return channel values
                if ~strcmp(varname,'ReachEleID') 
                    msg = sprintf('Reach property %s not available',varname);
                    warndlg(msg);
                    prop = []; return;
                end
            end
            %get channel ids from DispersionGraph based on channel connections
%             h_landward = Reach.LandwardPath(mobj); 
            estobj = getClassObj(mobj,'Inputs','Estuary');
            h_landward = type_sub_graph(estobj.DispersionGraph,'Channel');
            reachChannelID = h_landward.Nodes.EleID;            
            nrpt = length(reachChannelID);
            %DispersionGraphID and ReachELeID are variable length arrays so
            %return them as cells and everything else as an array
            if strcmp(varname,'DispersionGraphID') || strcmp(varname,'ReachEleID')
                prop = cell(nrpt,1);
                for i=1:nrpt
                    prop{i,1} = obj(i).(varname);
                end
            else
                prop = zeros(nrpt,1);
                for i=1:nrpt
                    prop(i,1) = obj(i).(varname);
                end
            end
        end
        
%%        
        function prop = getReachEleProp(mobj,varname)
            %property varname returned as an element array
            %get reach properties for each reach
            reachprop = Reach.getReachProp(mobj,varname);
            if isempty(reachprop), return; end
            %get ids of elements linked to channel in each reach
            reachids = Reach.getReachProp(mobj,'ReachEleID');
            %get channel ids from DispersionGraph based on channel connections
%             h_landward = Reach.LandwardPath(mobj);  
            estobj = getClassObj(mobj,'Inputs','Estuary');
            h_landward = type_sub_graph(estobj.DispersionGraph,'Channel');
            reachChannelID = h_landward.Nodes.EleID;
            %now assign reach property to each element
            nreach = length(reachChannelID);
            eleobj  = getClassObj(mobj,'Inputs','Element');
            prop = zeros(length(eleobj),1);
            for i=1:nreach
                prop(reachChannelID(i),1) = reachprop(i);
                if ~isempty(reachids)                    
                    idx = reachids{i};     %cell of ids for each reach
                    prop(idx,1) = reachprop(i);                    
                end
            end
            %deltas are outside and so not in a reach
            %assign them the properties of the linked channel element            
            eletype = getEleProp(eleobj,'EleType');
            extypes = mobj.GeoType(mobj.EXtypes);            
            idl = find(matches(eletype,extypes));
%             idl = find(contains(eletype,'Delta'));
            if ~isempty(idl)
                eleid = getEleProp(eleobj,'EleID');
                eleid = eleid(idl);
%                 estobj  = getClassObj(mobj,'Inputs','Estuary');
                h_reach = estobj.DispersionGraph;
                reachEleID = h_reach.Nodes.EleID;                
                for k=1:length(idl)
                    deltaidx= find(reachEleID==eleid(k)); %*idx ref h_reach ids
                    channelidx = successors(h_reach,deltaidx);
                    channelpidx = predecessors(h_reach,deltaidx);
                    if isempty(channelidx)
                        channelidx = channelpidx;
                    end
                    channelID = reachEleID(channelidx);   %*ID ref element id
                    idr = reachChannelID==channelID;      %idr ref h_landward ids
                    prop(idl(k),1) = reachprop(idr);
                end
            end
        end
        
%%
%         function h_landward = LandwardPath(mobj)
%             %the landward order of reaches is given by the channel subgraph
%             %of DispersionGraph 
%             estobj  = getClassObj(mobj,'Inputs','Estuary');
%             h_rch = estobj.DispersionGraph;                       
%             eletype = h_rch.Nodes.Type;
%             ich = find(strcmp(eletype,'Channel'));
%             h_landward = subgraph(h_rch,ich);            
%         end
        
%%        
%         function h_seaward = SeawardPath(mobj)
%             %seaward order of reaches is given by the inverse of the
%             %adjacency matrix for the landward order             
%             h_landward = Reach.LandwardPath(mobj);   
%             h_seaward = inverse_graph(h_landward);
%             nn = numnodes(h_landward);
%             [s,t] = findedge(h_landward);
%             pthAdj = sparse(s,t,h_landward.Edges.Weight,nn,nn);
%             h_seaward = digraph(pthAdj');
%             h_seaward.Nodes.EleID = h_landward.Nodes.EleID;
%             h_seaward.Nodes.Type  = h_landward.Nodes.Type;
%             h_seaward.Nodes.Name  = h_landward.Nodes.Name;
%         end
        
%%
%         function h_landward = LandwardFlow(h_flowpath)
%             %landward order of flow is given by the inverse of the
%             %adjacency matrix for the flowpath order  
%             nn = numnodes(h_flowpath);
%             [s,t] = findedge(h_flowpath);
%             pthAdj = sparse(s,t,h_flowpath.Edges.Weight,nn,nn);
%             h_landward = digraph(pthAdj');
%             h_landward.Nodes.EleID = h_flowpath.Nodes.EleID;
%             h_landward.Nodes.Type  = h_flowpath.Nodes.Type;
%             h_landward.Nodes.Name  = h_flowpath.Nodes.Name;
%         end
        
%%
        function Var = getReachEleVar(h_landward,EleVar)
            %assign elements variables to reaches
            rchID = h_landward.Nodes.EleID;
            eletype = h_landward.Nodes.Type;
            reachChannelID = rchID(strcmp(eletype,'Channel'));
            Var = EleVar(reachChannelID,:);       
        end
    
%%
        function AvVar = getAvUpstreamReachVar(h_landward,ReachVar,Eflag)
            %set path dependent parameters (uses graph ids)
            %average of reach and all adjacent upstream reach values
%             h_landward = Reach.LandwardFlow(h_flowgraph);  
%             h_landward = inverse_graph(h_flowgraph);
%             h_landward = type_sub_graph(h_landward,'Channel');
%             eletype = h_landward.Nodes.Type;
%             ich = find(strcmp(eletype,'Channel'));
%             h_landward = subgraph(h_landward,ich);             
            reachChannelID = h_landward.Nodes.EleID;
            nreach = length(reachChannelID);
            if Eflag    %input variable elements, so change to reaches
                ReachVar = Reach.getReachEleVar(h_landward,ReachVar);
            end
            AvVar = zeros(nreach,1);
            for i=1:nreach
                upstreamids = successors(h_landward,i);
                upVar = 0;
                if ~isempty(upstreamids)
                    for k=1:length(upstreamids)
                        upVar = upVar+ReachVar(upstreamids(k));
                    end
                end
                AvVar(i,1) = (ReachVar(i)+upVar)/2;
            end
            AvVar = AvVar(reachChannelID>0);
        end
        
%%
        function CumVar = getCumUpstreamReachVar(h_landward,ReachVar,Eflag)
            %set path dependent parameters (uses graph ids)
            %cumulative value of variable in upstream direction
%             h_landward = Reach.LandwardFlow(h_flowgraph);    
%             h_landward = inverse_graph(h_flowgraph);
            eletype = h_landward.Nodes.Type;
            ich = find(strcmp(eletype,'Channel'));
            h_landward = subgraph(h_landward,ich);             
            reachChannelID = h_landward.Nodes.EleID;
            nreach = length(reachChannelID);
            if Eflag    %input variable elements, so change to reaches
                ReachVar = Reach.getReachEleVar(h_landward,ReachVar);
            end
            CumVar = zeros(nreach,1);
            for i=1:nreach
                downstreamids = predecessors(h_landward,i);
                downVar = ReachVar(i);
                if ~isempty(downstreamids)
                    for k=1:length(downstreamids)
                        downVar = downVar+...
                            CumVar(downstreamids(k));
                    end
                end
                CumVar(i,1) = downVar;
            end
        end
        
%%
        function [network,answer] = getNetwork(caseDef,answer)
            %get the network to be used using either diffusion or flow
            % caseDef is an instance of mobj.Inputs or cobj.RunProps
            if nargin<2
                answer = questdlg('Select graph type','Network graph',...
                                'Network','Reach','Network');
            end
            disp = caseDef.Estuary.Dispersion;%dispersin matrix, d (m/s) - uses upper triangle 
            extdisp = caseDef.Estuary.ExternalDisp;
            eleobj = caseDef.Element;
            nodetxt = Estuary.setGraphVariables(eleobj);
            network = Estuary.setDispersionGraph(disp,extdisp,nodetxt); 

            switch answer
                case 'Network' %remove outside node
                    idc = find(~strcmp(network.Nodes.Type,'Outside'));             
                case 'Reach'    %select only channel nodes
                    idc = find(strcmp(network.Nodes.Type,'Channel'));
            end  
            network = subgraph(network,idc); 
        end   
%%
        function [reachlength,network] = getReachLengths(caseDef,ismidpt)
            %extract the distance between elements to construct along
            %channel plots as a function of distance
            % caseDef is an instance of mobj.Inputs or cobj.RunProps
            % ismidpt=true - lengths at mid-point of each element, otherwise
            % cumulative length to head of each element
            eLe = [caseDef.Element(:).Length]'; %Element.getEleProp(caseDef,'Length');
            network = Reach.getNetwork(caseDef,'Reach');
            CumX = Reach.getCumUpstreamReachVar(network,eLe,1);
            rchX = Reach.getReachEleVar(network,eLe);
            if ismidpt
                reachlength = CumX-rchX/2;  %take x as the mid-point of the reach            
            else
                reachlength = CumX;
            end
        end   
    end
%% ------------------------------------------------------------------------
% functions called by setReachProps
%--------------------------------------------------------------------------
    methods (Access=private)
        function setReachLength(obj,reachChannelID,h_landward)
            %set the reach length of each element and cumulative length from mouth
            for i=1:length(reachChannelID)
                downstreamids = predecessors(h_landward,i);
                cumlength = obj(i).ReachLength;
                if ~isempty(downstreamids)
                    for k=1:length(downstreamids)
                        cumlength = cumlength+obj(downstreamids(k)).CumulativeLength;
                    end
                end
                obj(i,1).CumulativeLength = cumlength;
            end
        end
%%
        
%%
        function [hw,mwl,lw,tr] = setHighLowWater(obj,reachChannelID,...
                                          rncobj,cstobj,advobj,wlvobj,damp)
            %set the high and low water levels using out put of CSTmodel if used
            if rncobj.IncDynHydraulics && isempty(cstobj)
                warndlg('Cannot include hydraulics because not defined')
                return
            elseif rncobj.IncDynHydraulics
                h_flowpath = advobj.RiverGraph;                 
                cstwls = assignCSTproperties(cstobj,obj,h_flowpath);
                %returns struct with mwl, amp, hwl and lwl
                cstmwl = cstwls.mwl;
            else
                cstmwl = 0;
            end
            hw =wlvobj.HWaterLevel*damp(reachChannelID);%reach specific HW
            lw = wlvobj.LWaterLevel*damp(reachChannelID);%reach specific LW 
            mwl = (wlvobj.MeanSeaLevel+cstmwl)*ones(size(reachChannelID));
            tr = (hw-lw);         %reach specific tidal range
        end
%%
        function setReachChannelProperties(obj,reachChannelID,reachEleID,...
                                                   eletype,h_landward)
            %set channel based properties of reach
            for i=1:length(reachChannelID) %uses graph ids and element ids
                %recursively get all elements linked to a reach channel
                idE = find(reachEleID==reachChannelID(i));
                nextElements = successors(h_landward,idE);
                nextElements = nextElements(~strcmp(eletype(nextElements),'Channel'));
                reachElements = nextElements;
                while ~isempty(nextElements)
                    nr = length(nextElements);
                    aele = [];
                    for jr=1:nr
                        sucs = successors(h_landward,nextElements(jr));
                        if ~isempty(sucs)
                            aele = [aele,sucs];  %#ok<AGROW>
                        end
                    end
                    nextElements = aele;
                    reachElements = [reachElements;nextElements];  %#ok<AGROW>
                end
                %linked elements for each reach based on DispersionGraph ids
                obj(i,1).DispersionGraphID = reachElements; 
                %linked elements for each reach based on Element ID, EleID
                idxReachEle = reachEleID(reachElements);
                obj(i,1).ReachEleID = idxReachEle;  
                obj(i,1).ReachChannelID = reachChannelID(i);
            end
        end
    %%
        function setCollectiveReachProperties(obj,reachChannelID,eleID,...
                                              eletype,V,S,n,HWL,MWL,LWL,TR)
            %%set collective properties for each reach. The properties set
            %include: HWvolume,HWarea,ReachPrism,ReachCSA,HWlevel,LWlevel
            for idr=1:length(reachChannelID) %uses graph ids and element ids
                idxReachEle = obj(idr,1).ReachEleID;
%                 reachElements = obj(i,1).DispersionGraphID;
                idc = find(eleID==reachChannelID(idr));
                obj(idr,1).HWvolume = getHWvolumes(obj,idr,idc,eletype,...
                                                                V,S,n,TR);
                
%                 if all(n(idxReachEle)>0)  %volumes are water volumes
%                     obj(idr,1).HWvolume = V(idc)+sum(V(idxReachEle))+S(idc)*TR(idr);
%                 else                      %volumes are sediment volumes
%                     idgt = n(idxReachEle)>0;
%                     idlt = 1-idgt;
%                     wetHWvols = sum(V(idxReachEle).*idgt);
%                     %tidal flat elements which are sediment volumes
%                     idfl = contains(eletype(reachElements),'flat',...
%                                                  'IgnoreCase',true).*idlt;
%                     sedHWflats = sum((S(idxReachEle)*TR(idr)-V(idxReachEle)).*idfl);
%                     %flood delta elements which are sediment volumes
%                     %this assumes that flood delta prism=S.tr, which
%                     %implies that the flood delta has minimal elevation
%                     %above low water (AsmitaOO used TR/2)
%                     idfd = contains(eletype(reachElements),'Flood').*idlt;
%                     sedHWfldelta = sum((S(idxReachEle)*TR(idr)).*idfd);
%                     sedHWvols = sedHWflats+sedHWfldelta;
%                     hwvol = V(idc)+wetHWvols+sedHWvols;          
%                     obj(idr,1).HWvolume = hwvol+S(idc)*TR(idr); %add prism volume over channel
%                 end   
                
                obj(idr,1).HWarea = S(idc)+sum(S(idxReachEle));
                obj(idr,1).ReachPrism = obj(idr).HWvolume-obj(idr).LWvolume;
                obj(idr,1).ReachCSA = (obj(idr).ReachPrism/2+ ...
                                     obj(idr).LWvolume)./obj(idr).ReachLength;   
                obj(idr,1).HWlevel = HWL(idr);
                obj(idr,1).MWlevel = MWL(idr); 
                obj(idr,1).LWlevel = LWL(idr);
            end
        end
%%
        function HWvolume = getHWvolumes(obj,idr,idc,eletype,V,S,n,TR)
            %determine the high water volume depending on use of water or
            %sediment volumes. idr is the reachChannelID and function is
            %called for each reach
            idxReachEle = obj(idr,1).ReachEleID;
            reachElements = obj(idr,1).DispersionGraphID;

            if all(n(idxReachEle)>0)  %volumes are water volumes
                hwvol = V(idc)+sum(V(idxReachEle)); %excludes volume over channel
            else                      %some volumes are sediment volumes
                idwet = n(idxReachEle)>0;  %index to water volumes
                idsed = ~idwet;            %index to sediment volumes 
                wetHWvols = sum(V(idxReachEle).*idwet);
                %tidal flat elements which are sediment volumes
                idfl = contains(eletype(reachElements),'flat',...
                                             'IgnoreCase',true).*idsed;
                sedHWflats = sum((S(idxReachEle)*TR(idr)-V(idxReachEle)).*idfl);
                %flood delta elements which are sediment volumes
                %this assumes that flood delta prism=S.tr, which
                %implies that the flood delta has minimal elevation
                %above low water (AsmitaOO used TR/2)
                idfd = contains(eletype(reachElements),'Flood').*idsed; %flood delta & sediment volume
                sedHWfldelta = sum((S(idxReachEle)*TR(idr)).*idfd);
                sedHWvols = sedHWflats+sedHWfldelta;
                hwvol = V(idc)+wetHWvols+sedHWvols;          
            end 
            HWvolume = hwvol+S(idc)*TR(idr); %add prism volume over channel
        end
%%
        function setUpstreamCSA(obj,reachChannelID,h_landward)
            %set path dependent upstream CSA (uses graph ids)
            for i=1:length(reachChannelID)
                upstreamids = successors(h_landward,i);
                upCSA = 0;
                if ~isempty(upstreamids)                    
                    for k=1:length(upstreamids)
                        upCSA = upCSA+obj(upstreamids(k)).ReachCSA;
                    end   
                end
                obj(i,1).UpstreamCSA = (obj(i).ReachCSA+upCSA)/2;
            end
        end
%%
        function setRiverFlow(obj,reachChannelID,advobj)
            %set river flow velocity at upstream end of each element
            riverDischarge = zeros(length(reachChannelID),1);
            if ~isempty(advobj) && ~isempty(advobj.RiverFlows) && any(advobj.RiverIn(:)~=0)
                %river flows have been defined and at least one input is defined
                h_flowpath = advobj.RiverGraph; 
                flowpathID = h_flowpath.Nodes.EleID;
                eletype = h_flowpath.Nodes.Type;
                ich = strcmp(eletype,'Channel');
                
                for i=1:length(flowpathID)
                    upstreamids = predecessors(h_flowpath,i);
                    if ~isempty(upstreamids)
                        idf = findedge(h_flowpath,upstreamids,i);
                        riverDischarge(i) = sum(h_flowpath.Edges.Weight(idf));
                    else
                        riverDischarge(i) = 9999; %river source (removed line 194)
                    end
                end
                riverDischarge = riverDischarge(ich);
            end   
            %now assign to reaches
            for idx=1:length(reachChannelID)
                obj(idx,1).RiverFlow = riverDischarge(idx)./obj(idx).UpstreamCSA;
            end 
        end
%%
        function setCumTidalPrism(obj,reachChannelID,h_landward)
            %cumulative upstream tidal prism (including reach prism of element)
            %if channels have been added out of order need to find upstream order 
            II = ones(length(reachChannelID),1);
            mouthids = find(logical(II-indegree(h_landward)));
            for j = 1:length(mouthids)
                reachOrder = dfsearch(h_landward,mouthids(j));            
                nreach = length(reachOrder);
                for i=1:nreach
                    idx = reachOrder(nreach+1-i);                
                    upstreamids = successors(h_landward,idx);
                    cumprism = obj(idx).ReachPrism;
                    if ~isempty(upstreamids)
                        for k=1:length(upstreamids)
                            cumprism=cumprism+obj(upstreamids(k)).UpstreamPrism;
                        end
                    end
                    obj(idx,1).UpstreamPrism = cumprism;
                end  
            end
        end        
    end
end