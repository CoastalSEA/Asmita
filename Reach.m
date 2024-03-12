classdef Reach < handle
    % Reaches are used to hold various gross properties that are used in the
    % model and results output. A reach is defined by any of the reach channel 
    % types, defined by RCtype (eg Channel), and any other reach elements,
    % defined by REtype, that are linked to the reach channel. RCtype and
    % REtype are properteis of the Asmita class that define subsets of GeoType
    %
    % NOTES
    % Graphs are used to define connectivity. 
    % * g_network defines the full connectivity as defined in the
    % DispersionGraph and includes an 'end node' for the Sea.
    % * g_landward is a subgraph of just the element types that define a 
    % reach (RCtypes) and has NO 'end nodes'.includes an 'end node' for the Sea.
    % * g_flowpath uses RiverGraph and includes all the elements in the 
    % advection path and includes 'end nodes' for the Sea and River(s).
    % Creating a subgraph using type_sub_graph(graph,ntype), or 
    % subgraph(graph,nodeids) are ways of removing these 'end nodes' to
    % enable direct indexing. Alternatively, the function ele2node finds
    % the graph node ids for given element ids or element type(s).
    % 
    % The Reach object array has an instance for each reach, including the
    % outside (end node).
    %
    %----------------------------------------------------------------------
    % AUTHOR
    % Ian Townend
    %
    % COPYRIGHT
    % CoastalSEA, (c) 2016
    %----------------------------------------------------------------------
    %  
    properties (Transient)
        ReachChannelID      %Element ID of Elements that define each reach
        ReachEleID          %Elements within a reach using Element ID
        HWvolume            %total volume of reach to high water
        LWvolume            %total volume to low water (usually just channel)
        HWarea              %total plan area of reach at high water
        LWarea              %total plan area of reach at low water
        HWlevel             %elevation of HW in reach at time t
        LWlevel             %elevation of LW in reach at time t
        MWlevel             %elevation of MWL in reach at time t
        dHWlevel            %change in elevation of HW in reach at time t
        dLWlevel            %change in elevation of LW in reach at time t
        dMWlevel            %change elevation of MWL in reach at time t
        TidalRange          %tidal range in reach (for use in output only)
        ReachPrism          %volume between high and low water in reach
        UpstreamPrism       %volume between high and low water of current 
                            %reach and all reaches upstream of current reach
        ReachLength         %length of reach (m)
        CumulativeLength    %distance to upstream end of element from mouth (excl.delta) (m)
        ReachCSA            %cross-sectional area of reach at mtl (m2)
        UpstreamCSA         %cross-sectional area of upstream end of reach at mtl (m2)
        RiverFlow           %Velocity at upstream end of reach based on river discharge and CSA
    end
    
    methods (Access=private)
        function obj = Reach
        end
    end
    
%%    
    methods (Static)
        function setReach(mobj,isprops)
            %initialise reach and assign ReachID to each element
            if nargin<2, isprops=false; end

            obj = Reach;                          %new instance
            %set ID properties of each reach (ReachEleID,ReachChannelID)
            obj = setReachChannelProperties(obj,mobj);
            
            %assign initial reach id to elements for use in plotting   
            eleobj = getClassObj(mobj,'Inputs','Element');
            for i=1:length(obj)
                %assign reach id to elements in landward order from Sea
                %elements that are not part of a reach are assigned id=0
                %this is based on element types that are not in REtypes.
                idc = obj(i).ReachChannelID;
                if idc>0
                    eleobj(idc).ReachID = idc;  %element is main reach element
                end
                idrele = obj(i).ReachEleID;
                for j=1:length(idrele)
                    eleobj(idrele(j)).ReachID = idc;
                end
            end            
            setClassObj(mobj,'Inputs','Element',eleobj);
            
            if isprops
                Reach.setReachProps(mobj); %also called directly in AsmitaModel.InitTimeStep
            end
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
            advobj = getClassObj(mobj,'Inputs','Advection');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');

            if isempty(eleobj)
                warndlg('No Elements defined');
                return
            elseif isempty(rncobj)
                warndlg('Run conditions not defined');
                return
            end
            
            if rncobj.IncRiver && ~isempty(advobj.RiverGraph)
                %use existing river graph
                g_flowpath = advobj.RiverGraph;
            elseif rncobj.IncRiver && ~isempty(advobj.RiverGraph) && ...
                ~isempty(advobj.RiverFlows) && any(advobj.RiverIn(:,2)~=0)
                %advection network has not been initialised
                g_flowpath = Advection.initialiseRiverGraph(mobj);
            else
                g_flowpath = [];
            end
            
            %set the LW channel based properties of reach
            obj = setLWproperties(obj,eleobj);

            [g_landward,~] = Reach.getReachGraph(mobj);
            %set reach length of each element and cumulative length from mouth
            setReachLength(obj,g_landward)

            %set high and low water levels using output of CSTmodel if used
            WL = setHighLowWater(obj,mobj,g_flowpath);

            %set collective properties for each reach (HWvolume,HWarea,
            %ReachPrism,ReachCSA,HWlevel,LWlevel,TidalRange)
            setCollectiveReachProperties(obj,mobj,WL);                                                    

            %set path dependent upstream CSA (uses graph ids)
            setUpstreamCSA(obj,g_landward);

            %river flow velocity at upstream end of element
            setRiverFlow(obj,mobj,g_flowpath);

            %cumulative upstream tidal prism (including reach prism of element)
            setCumTidalPrism(obj,g_landward)

            setClassObj(mobj,'Inputs','Reach',obj);
        end

%% ------------------------------------------------------------------------
% functions called to access properties from external methods
%--------------------------------------------------------------------------      
        function [g_landward,g_network] = getReachGraph(mobj)
            %returns network and reach graphs
            %uses the DispersionGraph to get the landward path, 
            estobj = getClassObj(mobj,'Inputs','Estuary');
            %graph of full network in landward direction (inc Sea)
            g_network = estobj.DispersionGraph; 
            if isempty(g_network)
                g_network = Estuary.initialiseDispersionGraph(mobj);                
            end
            rchtype = mobj.GeoType(mobj.RCtypes); %types that define a reach
            %graph of reaches in landward direction (inc Sea)
            g_landward = type_sub_graph(g_network,[{''},rchtype]); 
            %check that 'end node' for outside is still connected
            outids = outedges(g_landward,1);
            if isempty(outids)
                g_landward = addedge(g_landward,1,2,1);
            end
%             figure; plot(g_landward,'EdgeLabel',g_landward.Edges.Weight);
        end
%%
        function prop = getReachProp(mobj,varname,incall)
            %property varname returned as a reach array
            if nargin<3
                incall = false;  %incall - return all reaches if true,
            end                  %         omit outer reach if false
            
            obj = getClassObj(mobj,'Inputs','Reach');
            if isempty(obj)
                Reach.setReach(mobj,true);
            end

            % ReachELeID contains variable length arrays so
            %return as cells and everything else as an array
            nrpt = length(obj);
            if strcmp(varname,'ReachEleID')               
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
            %
            if ~incall
                prop = prop(2:end);
            end
        end
        
%%        
        function prop = getReachEleProp(mobj,varname)
            %property varname returned as an Element array, assigning reach
            %properties to Elements
            %get reach properties for each reach
            prop = [];
            reachprop = Reach.getReachProp(mobj,varname,true);
            if isempty(reachprop), return; end
            if strcmp(varname,'ReachEleID')
                warndlg('Use getReachProp to access ReachEleID')
                return;
            end
            %get ids of elements linked to channel in each reach
            reachids = Reach.getReachProp(mobj,'ReachEleID',true);
            reachChannelID = Reach.getReachProp(mobj,'ReachChannelID',true);
            nreach = length(reachChannelID);
            
            %now assign reach property to each element
            eleobj  = getClassObj(mobj,'Inputs','Element');
            prop = zeros(length(eleobj),1);
            for i=1:nreach
                if reachChannelID(i)>0
                    prop(reachChannelID(i),1) = reachprop(i);
                end
                idx = reachids{i};   %cell of ids for each reach
                if ~isempty(idx)                        
                    prop(idx,1) = reachprop(i);                    
                end
            end
        end
        
%%
        function Var = getReachEleVar(g_landward,EleVar)
            %assign Element properties, EleVar, to reaches
            idx = g_landward.Nodes.EleID>0;
            reachChannelID = g_landward.Nodes.EleID(idx);    
            if isvector(EleVar)
                Var = EleVar(reachChannelID);    %single property for each element
            else
                Var = EleVar(:,reachChannelID);  %rows are time cols are elements
            end
        end
    
%%
        function AvVar = getAvUpstreamReachVar(g_landward,ReachVar,Eflag)
            %set path dependent parameters (uses graph ids). average of 
            %ReachVar for each reach and all adjacent upstream reach values 
            %ReachVar can be an Element array if Eflag is true    
            if nargin<3, Eflag = false; end %assumes ReachVar input
            
            idx = g_landward.Nodes.EleID>0;
            reachChannelID = g_landward.Nodes.EleID(idx);  
            nreach = length(reachChannelID);
            if Eflag    %input variable elements, so change to reaches
                ReachVar = Reach.getReachEleVar(g_landward,ReachVar);
            end
            AvVar = zeros(nreach,1);
            for i=1:nreach
                upstreamids = successors(g_landward,i);
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
        function CumVar = getCumUpstreamReachVar(g_landward,ReachVar,Eflag)
            %get path dependent parameters (uses graph ids) cumulative value 
            %of ReachVar variable in upstream direction. ReachVar can be an 
            %Element array if Eflag is true    
            if nargin<3, Eflag = false; end %assumes ReachVar input
            
            idx = g_landward.Nodes.EleID>0;
            reachChannelID = g_landward.Nodes.EleID(idx);  
            nreach = length(reachChannelID);
            if Eflag    %input variable elements, so change to reaches
                ReachVar = Reach.getReachEleVar(g_landward,ReachVar);
            end
            CumVar = zeros(nreach,1);
            for i=1:nreach
                downstreamids = predecessors(g_landward,i);
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
        function [network,answer] = getNetwork(caseDef,answer)    %NOT TESTED
            %get the network to be used using either diffusion or flow
            % caseDef is an instance of mobj, mobj.Inputs or cobj.RunProps
            if nargin<2
                answer = questdlg('Select graph type','Network graph',...
                                'Network','Reach','Network');
            end
            
            inoutxt = {'Sea';'Rivers'};
            if isa(caseDef,'Asmita')
                eleobj = getClassObj(caseDef,'Inputs','Element');
                estobj = getClassObj(caseDef,'Inputs','Estuary');
                disp = estobj.Dispersion;       %dispersion matrix
                extdisp = estobj.ExternalDisp;  %external exchange
            else
                eleobj = caseDef.Element;
                disp = caseDef.Estuary.Dispersion;       %dispersion matrix
                extdisp = caseDef.Estuary.ExternalDisp;  %external exchange
            end
            nodetxt = setnodetext(eleobj,inoutxt);
            dispGraph = matrix2graph(disp,extdisp,[],nodetxt);

            switch answer
                case 'Network' %remove outside node
                    idc = find(~strcmp(dispGraph.Nodes.Type,''));             
                case 'Reach'   %select only reach nodes
                    idr = unique(getEleProp(eleobj,'ReachID'));
                    idc = ismember(dispGraph.Nodes.EleID,idr(idr>0));
            end  
            network = subgraph(dispGraph,idc); 
        end   
%%
        function [reachlength,network] = getReachLengths(caseDef,ismidpt)
            %extract the distance between elements to construct along
            %channel plots as a function of distance
            % caseDef is an instance of mobj.Inputs or cobj.RunParam
            % ismidpt=true - lengths at mid-point of each element, otherwise
            % cumulative length to head of each element
            network = Reach.getNetwork(caseDef,'Reach');
            eLe = [caseDef.Element(:).Length]';
            CumX = Reach.getCumUpstreamReachVar(network,eLe,true);
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
        function obj = setReachChannelProperties(obj,mobj)
            %set channel based properties of reach (ReachEleID,ReachChannel)

            %get network and landward graphs and asssociated IDs
            [g_landward,g_network]= Reach.getReachGraph(mobj);
            ID.networkEle = g_network.Nodes.EleID;    %Element IDs of the network graph
            ID.reachChannel = g_landward.Nodes.EleID; %Element IDs of the landward graph
            
            for i=1:length(ID.reachChannel) %uses graph ids and element ids
                %recursively get all elements linked to a reach channel
                idE = find(ID.networkEle==ID.reachChannel(i));
                nextElements = successors(g_network,idE);          
                nextReach = successors(g_landward,i); %connecting reach
                %remove elements that do not define a reach (eg EbbDelta)
                nextReach = validReachElements(obj,mobj,g_network,...
                                                 nextReach,'RCtypes');
                
                %handle external reach and check that reach elements are
                %valid members of set
                nextElements = validMemebers(obj,mobj,g_network,...
                                              nextElements,nextReach,ID,i);

                %recursive loop to find elements that are not directly
                %connected to the reach element (eg Saltmarsh via Tidalflat).
                reachElements = nextElements;  %linked elements in reach
                reachElements = linkedMembers(obj,mobj,g_network,...
                                              nextElements,reachElements);

                %linked elements for each reach based on Element ID, EleID
                idxReachEle = ID.networkEle(reachElements);
                obj(i,1).ReachEleID = idxReachEle;  
                obj(i,1).ReachChannelID = ID.reachChannel(i);
            end
            %returned as obj for use in setReach and saved to mobj.Inputs
            %for use in setReachProps
            setClassObj(mobj,'Inputs','Reach',obj);
        end
            
%%
        function nextElements = validMemebers(obj,mobj,g_network,...
                                            nextElements,nextReach,ID,ii)
            %handle external reach and check that reach elements are
            %valid members of set
            if ID.reachChannel(ii)==0 && ~isempty(nextReach)
                %for the 'end node' link any elements that do not form
                %a reach as defined by REtypes property in Asmita
                %NB this finds the nextElements that ARE reach elements
                %and uses the negation of this set in idx
                for i=1:length(nextReach)
                    idnextE = find(ID.networkEle==ID.reachChannel(nextReach(i)));
                    nextElements = [nextElements;successors(g_network,idnextE)];  %#ok<AGROW>

                    %remove elements that are not valid reach members
                    %including the reach channel element
                    nextChElements = validReachElements(obj,mobj,g_network,...
                                                 nextElements,'RCtypes');
                    nextRchElements = validReachElements(obj,mobj,g_network,...
                                                 nextElements,'REtypes');
                    nextReach = [nextReach;nextChElements;nextRchElements];   %#ok<AGROW>
                end
            elseif ID.reachChannel(ii)==0
                %the reach is an external element such as an ebb delta
            else
                %remove elements that are not valid reach members
                %including the reach channel element
                nextElements = validReachElements(obj,mobj,g_network,...
                                                nextElements,'REtypes');
            end                
            idx = ~ismember(nextElements,nextReach); %elements NOT in next reach
            nextElements = nextElements(idx);        %connecting elements
        end

%%
        function rchEles = validReachElements(~,mobj,g_network,...
                                                    linkedElements,type)
            %check that reach elements are valid members of the reach set
            %as defined by REtypes. 
            %type is RCtypes
            %linkedtypes is a list of element types that can define a reach
            %if type is REtypes
            %linkedtypes is a list of element types that are linked to a reach
            if strcmp(type,'RCtypes')
                retype = mobj.GeoType(mobj.RCtypes); %element types that define a reach
            else
                retype = mobj.GeoType(mobj.REtypes); %element types that link to a reach
            end
            linkedtypes = g_network.Nodes.Type(linkedElements);            
            ide = ismatch(linkedtypes,retype);   %elements that are reach elements
            rchEles = linkedElements(ide);
        end
%%
        function reachElements = linkedMembers(obj,mobj,g_network,...
                                                nextElements,reachElements)
             %recursive loop to find elements that are not directly
            %connected to the reach element (eg Saltmarsh via Tidalflat).
            while ~isempty(nextElements)
                nr = length(nextElements);
                aele = [];
                for jr=1:nr
                    sucs = successors(g_network,nextElements(jr));
                    if ~isempty(sucs)
                        aele = [aele,sucs];  %#ok<AGROW>
                    end
                end
                %remove elements that are not valid reach members
                %including the reach channel element
                nextElements = validReachElements(obj,mobj,g_network,...
                                                       aele,'REtypes');
                reachElements = [reachElements;nextElements];  %#ok<AGROW>
            end
        end
%%
        function obj = setLWproperties(obj,eleobj)
            %set the LW channel based properties of reach
            V = getEleProp(eleobj,'MovingVolume');
            S = getEleProp(eleobj,'SurfaceArea');
            eLe = getEleProp(eleobj,'Length');
            eleID = getEleProp(eleobj,'EleID');

            rchChID = [obj(:).ReachChannelID];
            for i=1:length(rchChID)                
                if rchChID(i)==0                    %set Sea reach to 0
                    obj(i,1).LWvolume = 0;      
                    obj(i,1).LWarea = 0;                
                    obj(i,1).ReachLength = 0;
                else                                %set reach values
                    idx = eleID==rchChID(i);
                    obj(i,1).LWvolume = V(idx);
                    obj(i,1).LWarea = S(idx);                
                    obj(i,1).ReachLength = eLe(idx);
                end
            end
        end
%%
        function setReachLength(obj,g_landward)
            %set the reach length of each element and cumulative length from mouth
            null = num2cell(zeros(1,length(obj)));  
            [obj.CumulativeLength] = null{:};  %reset cumulative length to zero
            rchChID = [obj(:).ReachChannelID]; %Element IDs of the reach channels          
           
            for i=1:length(rchChID)
                downstreamids = predecessors(g_landward,i);
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
        function WL = setHighLowWater(obj,mobj,g_flowpath)             
            %set the high and low water levels using output of CSTmodel if used
            eleobj = getClassObj(mobj,'Inputs','Element');
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            cstobj = getClassObj(mobj,'Inputs','CSThydraulics');

            rchChID = [obj(:).ReachChannelID]; %Element IDs of the reach channels
            damp = getEleProp(eleobj,'TidalDamping');
            unity = ones(size(rchChID));
            
            %initialise output arrays
            WL.hwl = unity*wlvobj.HWaterLevel;  %create defaults based on 
            WL.lwl = unity*wlvobj.LWaterLevel;  %values at the mouth
            WL.mwl = unity*wlvobj.MeanSeaLevel; %sets values for reach 0
            
            %adjust if dynamic hydraulics using CSTmodel is included
            rchChID = rchChID(2:end);
            if rncobj.IncDynHydraulics && isempty(cstobj)
                warndlg('Cannot include hydraulics because not defined')
                return
            elseif rncobj.IncDynHydraulics
                %returns struct with mwl, amp, hwl and lwl (mwl includes MTLatMouth)
                %adjust model results for any change in tidal amp and msl
                afact = (wlvobj.HWaterLevel-wlvobj.MeanSeaLevel)/wlvobj.TidalAmp;
                offset = wlvobj.MeanSeaLevel-wlvobj.MSL0;
                Q = g_flowpath.Edges.Weight(1);   %assumes all the same - ie single river input
                WL = assignCSTproperties(cstobj,obj,afact,offset,Q);
            else
                %assign reach values if no hydraulics
                ampdamp = (wlvobj.HWaterLevel-wlvobj.MeanSeaLevel)*damp(rchChID);
                WL.hwl(2:end) = wlvobj.MeanSeaLevel+ampdamp; %reach specific HW at time t
                WL.lwl(2:end) = wlvobj.MeanSeaLevel-ampdamp; %reach specific LW at time t
                WL.mwl(2:end) = wlvobj.MeanSeaLevel.*ones(size(rchChID));
            end
            WL.tr = (WL.hwl-WL.lwl);  %reach specific tidal range at time t
        end
    %%
        function setCollectiveReachProperties(obj,mobj,WL)
            %set collective properties for each reach. The properties set
            %include: HWvolume,HWarea,ReachPrism,ReachCSA,HWlevel,LWlevel,
            %TidalRange. uses graph ids (reachChannelID) and element ids (eleID)
            eleobj = getClassObj(mobj,'Inputs','Element');
            eleID = getEleProp(eleobj,'EleID');
            
            rchChID = [obj(:).ReachChannelID]; %Element IDs of the reach channels
            for idr=1:length(rchChID) 
                idc = find(eleID==rchChID(idr));
                if idc>0                   %reach area/volume properties
                    lobj = obj(idr);       %reach instance for idr
                    HWprops = getHWvolumes(lobj,mobj,idc,WL.tr(idr));                                                              
                    obj(idr,1).HWvolume = HWprops.volume;
                    obj(idr,1).HWarea = HWprops.area;
                    obj(idr,1).ReachPrism = lobj.HWvolume-lobj.LWvolume;
                    %reach average CSA at mtl
                    obj(idr,1).ReachCSA = (lobj.ReachPrism/2 + ...
                                        lobj.LWvolume)./lobj.ReachLength;                                       
                else                       %null values for external reach
                    obj(idr,1).HWvolume = 0;
                    obj(idr,1).ReachPrism = 0;
                    obj(idr,1).ReachCSA = 0;
                    obj(idr,1).HWarea = 0;
                end
                
                %assign water level properties to each reach
                obj(idr,1).dHWlevel = WL.hwl(idr)-obj(idr,1).HWlevel;
                obj(idr,1).dMWlevel = WL.mwl(idr)-obj(idr,1).MWlevel;
                obj(idr,1).dLWlevel = WL.lwl(idr)-obj(idr,1).LWlevel;
                obj(idr,1).HWlevel = WL.hwl(idr);
                obj(idr,1).MWlevel = WL.mwl(idr); 
                obj(idr,1).LWlevel = WL.lwl(idr);
                obj(idr,1).TidalRange =  WL.tr(idr);
            end
        end
%%
        function HWprops = getHWvolumes(obj,mobj,idc,TR)
            %determine the high water volume depending on use of water or
            %sediment volumes. idr is the reachChannelID and function is
            %called for each reach
            eleobj = getClassObj(mobj,'Inputs','Element');
            V = getEleProp(eleobj,'MovingVolume');
            S = getEleProp(eleobj,'SurfaceArea');
            n = getEleProp(eleobj,'TransportCoeff');
            eletype = getEleProp(eleobj,'EleType');
            idxReachEle = obj.ReachEleID;
            
            if all(n(idxReachEle)>0)  %volumes are water volumes
                hwvol = V(idc)+sum(V(idxReachEle)); %excludes volume over channel
            else                      %some volumes are sediment volumes
                idwet = n(idxReachEle)>0;  %index to water volumes
                idsed = ~idwet;            %index to sediment volumes 
                wetHWvols = sum(V(idxReachEle).*idwet);
                %tidal flat elements which are sediment volumes
                rchtypes = mobj.GeoType(mobj.REtypes); %reach element types
                flatypes =rchtypes(~strcmp(rchtypes,'FloodDelta')); %exclude FloodDelta
                idfl = ismatch(eletype(idxReachEle),flatypes).*idsed;
                sedHWflats = sum((S(idxReachEle)*TR-V(idxReachEle)).*idfl);
                %flood delta elements which are sediment volumes
                %this assumes that flood delta prism=S.tr, which
                %implies that the flood delta has minimal elevation
                %above low water (AsmitaOO used TR/2)
                idfd = ismatch(eletype(idxReachEle),'FloodDelta').*idsed; %flood delta & sediment volume
                sedHWfldelta = sum((S(idxReachEle)*TR).*idfd);
                sedHWvols = sedHWflats+sedHWfldelta;
                hwvol = V(idc)+wetHWvols+sedHWvols;          
            end 
            HWprops.volume = hwvol+S(idc)*TR; %add prism volume over channel
            HWprops.area = S(idc)+sum(S(idxReachEle));
        end
%%
        function setUpstreamCSA(obj,g_landward)
            %set path dependent upstream CSA at mtl (uses graph ids)
            %NB: this uses the reach average CSAs so that the model runs 
            %without the need to define the estuary convergence length, 
            %which is used for the tidal pumping estimates of CSA.
            rchChID = [obj(:).ReachChannelID]; %Element IDs of the reach channels
            for i=1:length(rchChID)
                upstreamids = successors(g_landward,i);
                upCSA = 0;
                if ~isempty(upstreamids)    
                    %sum the CSA of all the elements connected upstream
                    for k=1:length(upstreamids)
                        upCSA = upCSA+obj(upstreamids(k)).ReachCSA;
                    end   
                end
                %handle end reach for external elements
                if rchChID(i)>0  %reach element
                    obj(i,1).UpstreamCSA = (obj(i).ReachCSA+upCSA)/2;
                else             %Sea/external reach - rchChID=0
                    obj(i,1).UpstreamCSA = upCSA; %assign CSA of first reach 
                end                               %element to Sea reach
            end
        end
%%
        function setRiverFlow(obj,mobj,g_flowpath)
            %set river flow velocity at upstream end of each element
            rchChID = [obj(:).ReachChannelID]; %Element IDs of the reach channels
            
            if ~isempty(g_flowpath)
                riverDischarge = zeros(length(rchChID)-1,1);
                %river flows have been defined and at least one input is defined
                flowpathID = g_flowpath.Nodes.EleID(1:end-1);
                eletype = g_flowpath.Nodes.Type(1:end-1); %remove Rivers node                
                rchtypes = mobj.GeoType(mobj.RCtypes);    %reach element types
                ich = ismatch(eletype,rchtypes); 
                ich(1) = true;                            %add outside
                for i=1:length(flowpathID)
                    upstreamids = predecessors(g_flowpath,i);            
                    if ~isempty(upstreamids)
                        idf = findedge(g_flowpath,upstreamids,i);
                        riverDischarge(i) = sum(g_flowpath.Edges.Weight(idf));
                    end   
                end
                riverDischarge = riverDischarge(ich);
                %now assign to reaches
                riverflow = num2cell(riverDischarge./[obj(:).UpstreamCSA]');
            else
                riverflow = num2cell(zeros(1,length(rchChID)));
            end   
            
            [obj.RiverFlow] = riverflow{:};
        end
%%
        function setCumTidalPrism(obj,g_landward)
            %cumulative upstream tidal prism (including reach prism of element)
            %if channels have been added out of order need to find upstream order 
            rchChID = [obj(:).ReachChannelID]; %Element IDs of the reach channels
            II = ones(length(rchChID),1);
            mouthids = find(logical(II-indegree(g_landward)));
            for j = 1:length(mouthids)
                reachOrder = dfsearch(g_landward,mouthids(j));            
                nreach = length(reachOrder);
                for i=1:nreach
                    idx = reachOrder(nreach+1-i);                
                    upstreamids = successors(g_landward,idx);
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