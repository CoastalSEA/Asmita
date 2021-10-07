classdef River < handle
    %class to define river sources
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
        ChannelID = []        %EleID of channel river flows into 
        RiverFlow = 0         %flow rate from river at t=0 (m3/s) 
        RiverConc = 0         %concentration imported by advection (-)
        RiverTSC              %river timeseries of discharge (m3/s) and 
                              %sediment load (kg/m3)
    end
    
    properties (Transient)  %properties that are time varying
        tsRiverFlow           %flow rate from river at time,t (m3/s) 
        tsRiverConc           %concentration of imported sediment at time,t (-)
    end
    
    methods (Access=private)
        function obj = River
        end
    end
%%    
    methods (Static)
        function setRiver(mobj)
            %define or update properties for a river input
            eleobj  = getClassObj(mobj,'Inputs','Element');
            if isempty(eleobj)
                warndlg('Channel elements must be defined first');
                return
            end
            %
            obj  = getClassObj(mobj,'Inputs','River');
            if isempty(obj)
                obj = River;
                idx = 1;
            else
                elenames = getEleProp(eleobj,'EleName');
                [obj,idx] = editRiver(obj,elenames,'Add');
                if idx==0, return; end
            end
            ros = mobj.Constants.SedimentDensity;
            
            [prompt,~,defaultvalues] = riverProperties(obj,mobj);
            defaultvalues = defaultvalues(:,idx+1);
            title = 'River Input Parameters';
            numlines = 1;
            %use updated properties to call inpudlg and return new values
            ok = 0;
            while ok==0
                useInp=inputdlg(prompt,title,numlines,defaultvalues);
                if isempty(useInp)
                    return; %user cancels
                elseif str2double(useInp(1))==0 || isempty(useInp(1))  
                    return; %user removes flow input
                else
                    usenum(1) = str2double(useInp{1});
                    usenum(2) = str2double(useInp{2});    
                    usenum(3) = str2double(useInp{3});
                end 
                %check that river is assigned to a channel that exists
                eletype = getEleProp(eleobj,'EleType');
                ok = strcmp(eletype{usenum(1)},'Channel');
                if ok==0
                    warndlg('Can only connect to a channel, which must exit');
                end
            end
            
            obj(idx).ChannelID = usenum(1);
            obj(idx).RiverFlow = usenum(2);
            obj(idx).RiverConc = usenum(3)/ros;
            setClassObj(mobj,'Inputs','River',obj);
            
            %update Estuary Advection sources
            advobj = getClassObj(mobj,'Inputs','Advection');
            if isempty(advobj)          
                advobj = Advection.getAdvection;
            end
            %if empty initialise array for number of elements
            if isempty(advobj.RiverFlows)
                nele = length(eletype);
                advobj.RiverIn = zeros(nele,1);
            end
            rele = usenum(1);
            advobj.RiverIn(rele) = usenum(2);
            setClassObj(mobj,'Inputs','Advection',obj);
        end
%%
        function obj = addRiver(mobj,id,flow)
            %function to add a river when editing advection matrix
            %id is the channel id into which the river flows
            %flow is the flow rate of the river
            obj  = getClassObj(mobj,'Inputs','River');
            if isempty(obj)
                obj = River;
                idx = 1;
            else                
                obj = inp.River;
                idx = length(obj)+1;
                obj(idx,1) = River;
            end
            obj(idx).ChannelID = id;
            obj(idx).RiverFlow = flow;
            obj(idx).RiverConc = 0;
            setClassObj(mobj,'Inputs','River',obj);
        end 
%%
        function obj = delRiver(mobj)
            %delete a river input source
            msgtxt = 'No River inputs to delete';
            obj  = getClassObj(mobj,'Inputs','River',msgtxt);           
            if ~isempty(obj)
                eleobj  = getClassObj(mobj,'Inputs','Element');
                elenames = getEleProp(eleobj,'EleName');
                [obj,idx] = editRiver(obj,elenames,'List');
                if idx==0, return; end
                obj(idx) = [];
                setClassObj(mobj,'Inputs','River',obj);
                %update advection if set
                advobj  = getClassObj(mobj,'Inputs','Advection');
                if ~isempty(advobj)
                    rele = obj(idx).ChannelID;                
                    advobj.RiverIn(rele) = 0;
                    setClassObj(mobj,'Inputs','Advection',obj);
                end
            end
        end        
%%
       function prop = getRiverProp(mobj,varname)
           %function to access river properties from other functions
           obj  = getClassObj(mobj,'Inputs','River');
           eleobj  = getClassObj(mobj,'Inputs','Element');
           prop = zeros(length(eleobj),1);
           if isempty(obj) || isempty(obj(1).(varname))
                %no river or no variable defined
                return;
            end
            
            nriv = length(obj); %???????WHY NOT ASSIGNED BY ELEMENT            
            for i=1:nriv
                prop(i,1) = obj(i).(varname);
            end            
       end
%%   
        function [rele,okEle,msg,flow] = getSourceProps(mobj)
            %get the properties of a river source for advection            
            obj  = getClassObj(mobj,'Inputs','River');
            eleobj  = getClassObj(mobj,'Inputs','Element');
            eletype = getEleProp(eleobj,'EleType');
            okEle = strcmp(eletype,'Channel');                              
            msg{1} = 'Can only connect to a channel, which must exit\nRemove input to element %d';
            msg{2} = 'River input added to element %d\nProperties need to be defined in Setup>Rivers';
            msg{3} = 'River inputs have been added and need to be defined in Setup>Rivers';
            rele = []; flow = 0;
            if ~isempty(obj)
                rele = [obj.ChannelID]';
                flow = [obj.RiverFlow]';
            end
        end     
%%
        function obj = setRiverTS(mobj)
            %prompt user for file name and setup river timeseries
            %each river object can have an independent timeseries that is
            %loaded from a file (one river at a time)
            msgtxt = 'No rivers defined. Use River inputs to define some sources';
            obj  = getClassObj(mobj,'Inputs','River',msgtxt);  
            if isempty(obj)
                return;
            end
            %select river object to attach timeseries data
            %idx is the index of the river object (not the element)
            eleobj  = getClassObj(mobj,'Inputs','Element');
            elenames = getEleProp(eleobj,'EleName');
            if length(obj)>1
                [obj,idx] = editRiver(obj,elenames,'List');
                if idx==0, return; end
            else
                idx = 1;
            end
            RiverName = elenames{obj(idx).ChannelID};
            %get the name of the file and read it
            data = readinputdata();
            if isempty(data), return; end  %user aborted or file not read
            %
            myDatetime = data{1};
            obj(idx).RiverTSC = tscollection(myDatetime);
            obj(idx).RiverTSC.Name = 'River input';
            varnames = {'RiverFlow','SedLoad'};
            varunits = {'m3/year','kg/m3'};
            offset = 1;  %number of variables in 'data' already unpacked
            for i=offset+1:length(data)
                var = data{i};      
                % convert null records (99, 999, 9999) to NaN
                var(floor(var)==99 | var==99 | var==999 | var==9999) = NaN;
                tsv = timeseries(var,myDatetime);
                tsv.Name = varnames{i-1};
                tsv.DataInfo.Unit = varunits{i-1};
                tsv.UserData = sprintf('Flow into %s',RiverName);
                tsv.QualityInfo.Code = [0 1];
                tsv.QualityInfo.Description = {'good' 'bad'};

                obj(idx).RiverTSC = addts(obj(idx).RiverTSC,tsv);
            end     
            setClassObj(mobj,'Inputs','River',obj);
        end
        
%%
        function [flow,EleID] = getRiverTSprop(mobj,iriv,tsyear)
            %function to access time dependent river property from other functions
            % idrv is the index of the drift object
            % tsyear is the model time for which a drift is required
            obj  = getClassObj(mobj,'Inputs','River',msgtxt);  
            if isempty(obj) %no river defined
                flow = 0; EleID = 0;
                return;
            elseif isempty(obj(iriv).RiverTSC) %no timeseries input
                EleID = obj(iriv).ChannelID;
                flow = obj(iriv).RiverFlow;
                return;
            end

            EleID = obj(iriv).ChannelID;
            %interpolate time series to get new input value
            warning('off','MATLAB:linearinter:noextrap');                            
            ts1 = resample(obj(iriv).RiverTSC,tsyear);
            warning('on','MATLAB:linearinter:noextrap');
            tsvarnames = gettimeseriesnames(ts1);
            flow = ts1.RiverFlow.Data;
            if length(tsvarnames)<2
                conc = obj(iriv).RiverConc;
            else
                conc = ts1.SedLoad.Data/mobj.Constants.SedimentDensity;
            end
            
            %if outside data range, set default value                            
            if isnan(flow) || isempty(flow) 
                flow = obj(iriv).RiverFlow;
                conc = obj(iriv).RiverConc;
            end
            
            %update transient variables that are used in Model
            obj(iriv).tsRiverFlow = flow;
            obj(iriv).tsRiverConc = conc;
            setClassObj(mobj,'Inputs','River',obj);
        end
        
%%
        function TSplot(mobj,src,~)
            %plot the available river input timeseries data on a tab 
            ht = findobj(src,'Type','Axes');
            delete(ht);
            
            msgtxt = 'No River inputs defined';
            obj  = getClassObj(mobj,'Inputs','River',msgtxt);           
            if ~isempty(obj)
                return;
            else
                tsid = zeros(length(obj),1);
                numtimeseries = 0;
                for i=1:length(obj)
                    tsid(i) = ~isempty(obj(i).RiverTSC);
                    numtimeseries = numtimeseries + tsid(i);                    
                end
                if  numtimeseries<1
                    warndlg('No river timeseries defined')
                    return;
                end
            end
            
            axes('Parent',src);
            sp1 = subplot(2,1,1);
            tsid = find(tsid>0);
            for i=1:numtimeseries
                eleid = obj(tsid(i)).ChannelID;
                eleobj  = getClassObj(mobj,'Inputs','Element');
                elename = eleobj(eleid).EleName;
                legtxt = sprintf('Input to %s',char(elename));
                ts1 = obj(tsid(i)).RiverTSC.RiverFlow;            
                hp(i) = plot(ts1,'DisplayName',legtxt); %#ok<AGROW>
                if numtimeseries>1 && i==1
                    hold on; 
                end                
            end
            hold off
            title('','Visible','off');
            xlabel('','Visible','off');
            ylabel('River discharge (m^3/s)');
            legend(sp1,hp,'Location','best');
            
            sp2 = subplot(2,1,2);
            for i=1:numtimeseries
                eleid = obj(tsid(i)).ChannelID;
                elename = inp.Element(eleid).EleName;
                legtxt = sprintf('Input to %s',char(elename));
                tsvarnames = gettimeseriesnames(obj(tsid(i)).RiverTSC);
                if length(tsvarnames)>1
                    ts1 = obj(tsid(i)).RiverTSC.SedLoad;            
                    gp(i) = plot(ts1,'DisplayName',legtxt); %#ok<AGROW>
                else
                    yval = obj(tsid(i)).RiverConc;
                    yval = yval*mobj.Constants.SedimentDensity;
                    gp(i) = plot([min(hp(1).XData),max(hp(1).XData)],...
                        [yval,yval],'DisplayName','Constant value'); %#ok<AGROW>
                end
                if numtimeseries>1 && i==1
                    hold on; 
                end                
            end
            hold off            
            title('','Visible','off');
            xlabel('Year');
            ylabel('Sediment load (kg/m^3)');
            legend(sp2,gp,'Location','best');
        end
        
    end    
%%    
    methods
        function initialiseFlow(obj)
            %set the transient values at the start of a model run
            for i=1:length(obj)
                obj(i).tsRiverFlow = obj(i).RiverFlow;       
                obj(i).tsRiverConc = obj(i).RiverConc;  
            end
        end
%%        
        function riverPropertiesTable(obj,src,mobj)
            %table for Input Summary tab display of river properties
            [~,colnames,userdata] = riverProperties(obj,mobj);
            numrivers = length(obj);
            cwidth = {160 45};
            for j =1:numrivers
                cwidth{j+1} = 45;  %add column widths based on no. sources
            end  
            %add row for inclusion of timeseries
            [nrow,~] = size(userdata);
            userdata{nrow+1,1} = 'Timeseries included';
            for i=1:length(obj)
                if isempty(obj(i).RiverTSC)
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
            Table.Position(2)=0.85-Table.Extent(4); 
        end
    end
    
%%
    methods (Access=private)
        function [obj,idx] = editRiver(obj,elenames,flag)
            %edit the properties of a river instance            
            numRivers = length(obj);
            riversList = cell(numRivers,1);
            for i=1:numRivers
                channelid = num2str(obj(i).ChannelID);
                riversList{i} = sprintf('%s: %s',channelid,elenames{obj(i).ChannelID});
            end
            if strcmp(flag,'Add')
                riversList{numRivers+1} = 'Add a River';
            end
            
            if ~isempty(riversList)
                [idx, ok] = listdlg('Name','River Inputs', ...
                    'PromptString','Select Channel for river input', ...
                    'SelectionMode','single', ...
                    'ListString',riversList);
                if ok<1, idx=0; end
            else
                idx = 0;
            end
            
            if idx>numRivers
                idx = numRivers+1;  %add a river
                obj(idx,1) = River;
            end
        end
        
%%        
        function [rownames,colnames,userdata] = riverProperties(obj,mobj)
            %define river Properties  
            ros = mobj.Constants.SedimentDensity;
            rownames = {'Channel ID for river input',...
                        'Flow rate (m^3/s)',...
                        'Sediment conc. density in river (kg/m3)'};
            colnames = {'Property','Values'};  
            numrivers = length(obj);
            userdata = cell(length(rownames),numrivers+1);
            for i=1:length(rownames)
                userdata{i,1} = rownames{i};
            end   
            for j = 1:numrivers
                k = j+1;
                userdata{1,k} = num2str(obj(j).ChannelID);
                userdata{2,k} = num2str(obj(j).RiverFlow);
                if isfield(mobj.Inputs.Advection,'RiverIn') && ...
                            ~isempty(mobj.Inputs.Advection.RiverIn)
                    userdata{2,k} = num2str(inp.Advection.RiverIn(obj(j).ChannelID));
                end
                rivConc = obj(j).RiverConc*ros;
                userdata{3,k} = num2str(rivConc);   
            end
        end
        
    end
end