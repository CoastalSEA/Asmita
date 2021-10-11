classdef Drift < handle
    % class to define littoral drift sources
    %
    %----------------------------------------------------------------------
    % AUTHOR
    % Ian Townend
    %
    % COPYRIGHT
    % CoastalSEA, (c) 2016
    %----------------------------------------------------------------------
    %
    properties
        DriftEleID = []       %EleID of element drift flows into  
        DriftRate = 0         %littoral drift rate as volume of sediment (m3/year)
        DriftTSC              %timeseries collection of drift rate (m3/year)
    end
    
    properties (Transient)  %properties that are time varying
        tsDriftRate           %drift rate at time,t (m3/s) 
    end

    methods (Access=private)
        function obj = Drift
        end
    end
%%    
    methods (Static)
        function setDrift(mobj)
            %define or update properties for a drift input
            eleobj  = getClassObj(mobj,'Inputs','Element');
            if isempty(eleobj)
                warndlg('Channel elements must be defined first');
                return
            end
            %
            obj  = getClassObj(mobj,'Inputs','Drift');
            if isempty(obj)
                obj = Drift;
                idx = 1;
            else
                [obj,idx] = getDriftSource(obj,eleobj,'Add');
                if idx==0, return; end
            end

            [prompt,~,defaultvalues] = driftProperties(obj,mobj);
            defaultvalues = defaultvalues(:,idx+1);
            title = 'Drift Input Parameters';
            numlines = 1;
            %use updated properties to call inpudlg and return new values
            ok = 0;
            while ok==0
                useInp=inputdlg(prompt,title,numlines,defaultvalues);
                if isempty(useInp)
                    return; %user cancels
                elseif str2double(useInp{1})==0 || isempty(useInp{1})
                    return; %user removes flow input
                else
                    usenum(1) = str2double(useInp{1});
                    usenum(2) = str2double(useInp{2});
                end 
                %check that drift is assigned to a channel that exists
                eletype = getEleProp(eleobj,'EleType');
                ok = contains(eletype{usenum(1)},'Delta');
                ok = ok+contains(eletype{usenum(1)},'Beach');
                ok = ok+contains(eletype{usenum(1)},'Shore');
                ok = ok+contains(eletype{usenum(1)},'Spit');
                if ok==0
                    warndlg('Can only connect to delta, beach, shore, or spit which must exit');
                end
            end            
            obj(idx).DriftEleID = usenum(1);
            obj(idx).DriftRate = usenum(2);
            setClassObj(mobj,'Inputs','Drift',obj);

            %update Estuary Advection sources
            advobj = getClassObj(mobj,'Inputs','Advection');
            if isempty(advobj)          
                advobj = Advection;
            end
            %if empty initialise array for number of elements
            if isempty(advobj.DriftFlows)
                nele = length(eletype);
                advobj.DriftIn = zeros(nele,1);
            end
            rele = [eleobj(:).EleID]==usenum(1);
            advobj.DriftIn(rele) = usenum(2);
            
            setClassObj(mobj,'Inputs','Advection',advobj); 
            warndlg({'Drift ADDED';
                     'Do not forget to update Drift Advection'});
        end
%%
        function addDrift(mobj,id,flow)
            %function to add a drift when editing Advection matrix
            %id is the channel id into which the drift flows
            %flow is the littoral drift rate
            obj  = getClassObj(mobj,'Inputs','Drift');
            if isempty(obj)
                obj = Drift;
                idx = 1;                
            else
                idx = length(obj)+1;
                obj(idx,1) = Drift;
            end
            obj(idx).DriftEleID = id;
            obj(idx).DriftRate = flow;
            
            setClassObj(mobj,'Inputs','Drift',obj);
            warndlg({'Drift ADDED';
                     'Do not forget to update Drift Advection'});
        end
%%
        function delDrift(mobj)
            %delete a drift input source
            msgtxt = 'No Drift inputs to delete';
            obj  = getClassObj(mobj,'Inputs','Drift',msgtxt);           
            if isempty(obj), return; end
            
            eleobj  = getClassObj(mobj,'Inputs','Element');
            [obj,idx] = getDriftSource(obj,eleobj,'List');
            if idx==0, return; end
            id_driftinput = obj(idx).DriftEleID;
            obj(idx) = [];
            setClassObj(mobj,'Inputs','Drift',obj);

            %update advection if set
            advobj  = getClassObj(mobj,'Inputs','Advection');
            if ~isempty(advobj)
                idele = [eleobj(:).EleID]==id_driftinput;
                advobj.DriftIn(idele) = 0;                    
                setClassObj(mobj,'Inputs','Advection',advobj);
            end
            
            warndlg({'Drift DELETED';
                     'Do not forget to update Drift Advection'});
        end        
%%
       function prop = getDriftProp(mobj,varname)
           %function to access drift properties from other functions
           obj  = getClassObj(mobj,'Inputs','Drift');
           prop = zeros(length(mobj.h_ele),1);
           if isempty(mobj.h_dft) || isempty(mobj.h_dft(1).(varname))
               %no drift of no variable defined
               return;
           end
           
           eleobj  = getClassObj(mobj,'Inputs','Element');
           prop = zeros(length(eleobj),1);
           ndft = length(obj);
           for i=1:ndft
               %assign value to the element the drift source flows into
               idele = [eleobj(:).EleID]==obj(i).DriftEleID;
               prop(idele,1) = obj(i).(varname);
           end
       end
%%   
        function [rele,okEle,msg,flow] = getSourceProps(mobj)
            %get the properties of a drift source for advection
            obj  = getClassObj(mobj,'Inputs','Drift');
            eleobj  = getClassObj(mobj,'Inputs','Element');
            eletype = getEleProp(eleobj,'EleType');
            okEle = strcmp(eletype,'Delta');
            okEle = okEle+strcmp(eletype,'Beach');
            okEle = okEle+strcmp(eletype,'Shoreface');
            msg{1} = 'Can only connect to delta/beach/shoreface elements, which must exit\nRemove input to element %d';
            msg{2} = 'Drift input added to element %d\nProperties need to be defined in Setup>Drift';
            msg{3} = 'Drift inputs have been added';            
            ndft = length(obj);
            rele(ndft,1) = 0; flow(ndft,1) = 0;
            for i=1:ndft
                rele(i) = find([eleobj(:).EleID]==obj(i).DriftEleID);
                flow(i) = obj(i).DriftRate;
            end
        end  
%%
        function setDriftTS(mobj)
            %prompt user for file name and setup drift timeseries
            %each drift object can have an independent timeseries that is
            %loaded from a file. This can define either just the input drift,
            %or drift between selected elements in the flow path
            msgtxt = 'No drifts defined. Use Drift inputs to define some sources';
            obj  = getClassObj(mobj,'Inputs','Drift',msgtxt);  
            if isempty(obj)
                return;
            end
            %select which drift object the file is to be attached to
            %idx is the index of the drift object (not the element)
            eleobj  = getClassObj(mobj,'Inputs','Element');
            elenames = getEleProp(eleobj,'EleName');
            if length(obj)>1                
                [obj,idx] = getDriftSource(iobj,elenames,'List');
                if idx==0, return; end
            else
                idx = 1;                
            end
            
            DriftName = elenames{obj(idx).DriftEleID};            
            %get the name of the file and read it
            [data,header] = readinputdata();
            if isempty(data), return; end  %user aborted or file not read
            %use header definitions to creat timeseries collection
            pathtxt = split(header);
            pathtxt = pathtxt(2:end,1); %remove date string
            pathSrt = zeros(length(pathtxt),1); pathEnd = pathSrt;
            for i=1:length(pathtxt)
                pathstr = char(pathtxt(i,1));
                startid = regexp(pathstr,'(')+1;
                endid = regexp(pathstr,')')-1;
                if ~isempty(startid)
                    temp = split(pathstr(startid:endid),',');
                    pathSrt(i,1) = str2double(temp(1)); %array of start nodes
                    pathEnd(i,1) = str2double(temp(2)); %array of end nodes
                end
            end
            %
            myDatetime = data{1};
            %
            obj(idx).DriftTSC = tscollection(myDatetime);
            obj(idx).DriftTSC.Name = 'Drift';
            offset = 1;  %number of variables in 'data' already unpacked
            for i=offset+1:length(data)
                var = data{i};      
                % convert null records (99, 999, 9999) to NaN
                var(floor(var)==99 | var==99 | var==999 | var==9999) = NaN;
                tsv = timeseries(var,myDatetime);
                tsv.Name = sprintf('Qs_%d_%d',pathSrt(i-1),pathEnd(i-1));
                tsv.DataInfo.Unit = 'm3/year';
                tsv.UserData = sprintf('Drift starting from %s',DriftName);
                tsv.QualityInfo.Code = [0 1];
                tsv.QualityInfo.Description = {'good' 'bad'};

                obj(idx).DriftTSC = addts(obj(idx).DriftTSC,tsv);
            end
            
           setClassObj(mobj,'Inputs','Drift',obj);
        end
%%
        function [prop,EleID] = getDriftTSprop(mobj,idft,tsyear)
            %function to access time dependent drift property from other functions
            % idft is the index of the drift object
            % tsyear is the model time for which a drift is required
            obj  = getClassObj(mobj,'Inputs','Drift',msgtxt);  
            if isempty(obj) %no drift defined
                prop = 0; EleID = 0;
                return;
            elseif isempty(obj(idft).DriftTSC) %no timeseries input
                EleID = obj(idft).ChannelID;
                prop = obj(idft).DriftRate;
                return;
            end
            EleID = obj(idft).DriftEleID;
            %interpolate time series to get new input value
            warning('off','MATLAB:linearinter:noextrap');                            
            ts1 = resample(obj(idft).DriftTSC,tsyear);
            warning('on','MATLAB:linearinter:noextrap'); 
            varnames = gettimeseriesnames(ts1);
            nvar = length(varnames);
            if nvar>1
                EleID = unpackVarIDs(varnames);
            end
            prop = zeros(nvar,1);
            for k=1:length(varnames)
                indrift = ts1.(varnames{k}).Data;
                %if outside data range, set default value  
                if isnan(indrift) || isempty(indrift)
                    prop(k) = obj(idft).DriftRate;
                else
                    prop(k) = indrift;
                end
            end
            %update transient variables that are used in Model
            obj(idft).tsDriftRate = prop;
            
            setClassObj(mobj,'Inputs','Drift',obj);
            %--------------------------------------------------------------
            function EleID = unpackVarIDs(varnames)
                %recover the node ids for each edge defined in tsc
                EleID = zeros(nvar,2);
                for i = 1:nvar
                    avariable = varnames{i};
                    eleids = split(avariable(4:end),'_');
                    EleID(i,:) = str2double(eleids)';
                end                
            end
        end
        
%%
        function TSplot(mobj,src,~)
            %plot the available drift input timeseries data on a tab  
            ht = findobj(src,'Type','Axes');
            delete(ht);
            
            msgtxt = 'No Drift inputs defined';
            obj  = getClassObj(mobj,'Inputs','Drift',msgtxt);           
            if ~isempty(obj)
                return;
            else
                tsid = zeros(length(obj),1);
                numtimeseries = 0;
                for i=1:length(obj)
                    tsid(i) = ~isempty(obj(i).DriftTSC);
                    numtimeseries = numtimeseries + tsid(i);                    
                end
                if  numtimeseries<1
                    warndlg('No drift timeseries defined')
                    return;
                end
            end

            figax = axes('Parent',src);
            tsid = find(tsid>0);
            k = 0;
            for i=1:numtimeseries
                eleid = obj(tsid(i)).DriftEleID;
                tsc1 = obj(tsid(i)).DriftTSC;
                tsnames = gettimeseriesnames(tsc1);
                for j = 1:length(tsnames)
                    k = k+1;                    
                    hp(k) = plot(tsc1.(tsnames{j}),'DisplayName',tsnames{j}); %#ok<AGROW>
                    if numtimeseries>1 || length(tsnames)>1 && k==1 
                        hold on; 
                    end   
                end 
            end
            hold off
            
            if length(tsnames)>1
                legtxt = 'Drift input to multiple elements';
            else
                eleobj  = getClassObj(mobj,'Inputs','Element');
                elename = eleobj(eleid).EleName;
                legtxt = sprintf('Drift input to %s',char(elename));
            end
            xlabel('Year');
            ylabel('Drift rate (m^3/year)');
            title(legtxt);
            legend(figax,hp,'Location','best');
        end
    end
%%    
    methods
        function initialiseFlow(obj)
            %set the transient values at the start of a model run
            for i=1:length(obj)
                obj(i).tsDriftRate = obj(i).DriftRate;        
            end
        end
%%        
        function driftPropertiesTable(obj,src,mobj)
            %table for Input Summary tab display of drift properties
            [~,colnames,userdata] = driftProperties(obj,mobj);
            numrdrifts = length(obj);
            cwidth = {160 45};
            for j =1:numrdrifts
                cwidth{j+1} = 45;  %add column widths based on no. sources
            end  
            %add row for inclusion of timeseries
            [nrow,~] = size(userdata);
            userdata{nrow+1,1} = 'Timeseries included';
            for i=1:length(obj)
                if isempty(obj(i).DriftTSC)
                    userdata{nrow+1,i+1} = 'No';
                else
                    userdata{nrow+1,i+1} = 'Yes';
                end
            end
            
            Table = uitable('Parent',src, ...
                'ColumnName', colnames, ...
                'RowName', [], ....
                'ColumnWidth', cwidth, ...
                'Data',userdata, ...
                'Units','normalized');
            Table.Position(3:4)=Table.Extent(3:4);
            Table.Position(2)=0.4-Table.Extent(4); 
        end
    end    

%%
    methods (Access=private)
        function [obj,idx] = getDriftSource(obj,eleobj,flag)
            %edit the properties of a drift instance
            numDrifts = length(obj);
            driftList = cell(numDrifts,1);
            elenames = getEleProp(eleobj,'EleName');
            for i=1:numDrifts
                drifteleid = num2str(obj(i).DriftEleID);
                idele = [eleobj(:).EleID]==obj(i).DriftEleID;
                driftList{i} = sprintf('%s: %s',drifteleid,elenames{idele});
            end
            if strcmp(flag,'Add')
                driftList{numDrifts+1} = 'Add a Drift';
            end
            
            if ~isempty(driftList)
                [idx, ok] = listdlg('Name','Drift Inputs', ...
                    'PromptString','Select Element for drift input', ...
                    'SelectionMode','single', ...
                    'ListString',driftList);
                if ok<1, idx=0; end
            else
                idx = 0;
            end
            
            if idx>numDrifts
                idx = numDrifts+1;  %add a drift
                obj(idx,1) = Drift;
            end
        end
        
%%        
        function [rownames,colnames,userdata] = driftProperties(obj,mobj)
            %define drift Properties 
            advobj = getClassObj(mobj,'Inputs','Advection');
            rownames = {'Element ID for drift input',...
                        'Drift rate (m^3/year)'};
            colnames = {'Property','Value'};  
            numdrifts = length(obj);
            userdata = cell(length(rownames),numdrifts+1);
            for i=1:length(rownames)
                userdata{i,1} = rownames{i};
            end   

            for j = 1:numdrifts
                k = j+1;
                userdata{1,k} = num2str(obj(j).DriftEleID);
                userdata{2,k} = num2str(obj(j).DriftRate);
                if isfield(advobj,'DriftIn') && ~isempty(advobj.DriftIn)                                    
                    userdata{2,k} = num2str(inp.Advection.DriftIn(obj(j).DriftEleID));
                end   
            end
        end                
        
    end
end