classdef River < matlab.mixin.Copyable
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
        RiverRho = 0          %concentration density imported by advection (kg/m3)
        RiverTSC              %river timeseries of discharge (m3/s) and 
                              %sediment load (kg/m3)
    end
    
    properties (Transient)    %properties that are time varying
        tsRiverFlow           %flow rate from river at time,t (m3/s) 
        tsRiverConc           %concentration of imported sediment at time,t (-)
    end
    
    properties (Dependent)
        RiverConc             %concentration imported by advection (-)
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
                [obj,idx] = editRiver(obj,eleobj,'Add');
                if idx==0, return; end
            end
             
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
            obj(idx).RiverRho = usenum(3);
            setClassObj(mobj,'Inputs','River',obj);
            
            %update Estuary Advection sources
            advobj = getClassObj(mobj,'Inputs','Advection');
            if isempty(advobj)          
                advobj = Advection;
            end
            %if empty initialise array for number of elements
            if isempty(advobj.RiverFlows)
                nele = length(eletype);
                advobj.RiverIn = zeros(nele,2);
            end
            rele = [eleobj(:).EleID]==usenum(1);
            advobj.RiverIn(rele,2) = usenum(2);
            
            setClassObj(mobj,'Inputs','Advection',advobj);
            getdialog({'River ADDED';
                     'Do not forget to update River Advection'});
        end
%%
        function addRiver(mobj,id,flow)
            %function to add a river when editing advection matrix
            %id is the channel id into which the river flows
            %flow is the flow rate of the river
            obj  = getClassObj(mobj,'Inputs','River');
            if isempty(obj)
                obj = River;
                idx = 1;
            else                
                idx = length(obj)+1;
                obj(idx,1) = River;
            end
            obj(idx).ChannelID = id;
            obj(idx).RiverFlow = flow;
            obj(idx).RiverConc = 0;
            
            setClassObj(mobj,'Inputs','River',obj);
            getdialog({'River ADDED';
                     'Do not forget to update River Advection'});
        end 
%%
        function delRiver(mobj)
            %delete a river input source
            msgtxt = 'No River inputs to delete';
            obj  = getClassObj(mobj,'Inputs','River',msgtxt);           
            if isempty(obj), return; end
            
            eleobj  = getClassObj(mobj,'Inputs','Element');                               
            [obj,idx] = editRiver(obj,eleobj,'List');
            if idx==0, return; end
            id_riverinput = obj(idx).ChannelID;
            obj(idx) = [];
            setClassObj(mobj,'Inputs','River',obj);

            %update advection if set
            advobj  = getClassObj(mobj,'Inputs','Advection');
            if ~isempty(advobj)
                idele = [eleobj(:).EleID]==id_riverinput;
                advobj.RiverIn(idele) = 0;
                setClassObj(mobj,'Inputs','Advection',advobj);
            end

            getdialog({'River DELETED';
                     'Do not forget to update River Advection'});
        end        
%%
       function prop = getRiverProp(mobj,varname)
           %function to access river properties from other functions
           eleobj  = getClassObj(mobj,'Inputs','Element');
           prop = zeros(length(eleobj),1);
           obj  = getClassObj(mobj,'Inputs','River');           
           if isempty(obj) || isempty(obj(1).(varname))
               %no river or no variable defined
               return;
           end

           nriv = length(obj);
           for i=1:nriv
               %assign value to the element the river source flows into
               idele = [eleobj(:).EleID]==obj(i).ChannelID;
               prop(idele,1) = obj(i).(varname);
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
            % idrv is the index of the river object
            % tsyear is the model time for which a river is required
            obj  = getClassObj(mobj,'Inputs','River');  
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
            obj  = getClassObj(mobj,'Inputs','River');           
            if isempty(obj)
                getdialog(msgtxt);
                return;
            else
                tsid = zeros(length(obj),1);
                numtimeseries = 0;
                for i=1:length(obj)
                    tsid(i) = ~isempty(obj(i).RiverTSC);
                    numtimeseries = numtimeseries + tsid(i);                    
                end
                if  numtimeseries<1
                    getdialog('No river timeseries defined')
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
                elename = eleobj(eleid).EleName;
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
    methods (Static,Hidden)
        function setNewRiver(mobj)
            %initialise an empty instance of River (used in asm_oo2mui)
            obj = River;
            setClassObj(mobj,'Inputs','River',obj);
        end
    end    
%%    
    methods
        function eqConc = get.RiverConc(obj)
            %dependent property derived from RiverRho
            cn = muiConstants.Evoke;
            eqConc = obj.RiverRho/cn.SedimentDensity;
        end  
%%
        function initialiseFlow(obj)
            %set the transient values at the start of a model run
            for i=1:length(obj)
                obj(i).tsRiverFlow = obj(i).RiverFlow;       
                obj(i).tsRiverConc = obj(i).RiverConc;  
            end
        end
%%   
        function source = getSourceProps(obj,eleobj,~)
            %get the properties of a river source for advection            
            eletype = getEleProp(eleobj,'EleType');
            source.okEle = strcmp(eletype,'Channel');                              
            msg{1} = 'Can only connect to a channel, which must exit\nRemove input to element %d';
            msg{2} = 'River input added to element %d\nProperties need to be defined in Setup>Rivers';
            msg{3} = 'River inputs have been added and need to be defined in Setup>Rivers';
            nriv = length(obj);
            rele(nriv,1) = 0; flow(nriv,1) = 0;
            for i=1:nriv
                rele(i) = find([eleobj(:).EleID]==obj(i).ChannelID);
                flow(i) = obj(i).RiverFlow;
            end
            source.InputEle = rele; 
            source.msg = msg; 
            source.flow = flow;
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
%%
        function tabPlot(obj,src,mobj)
            %plot tidal pumping advection tab
            msgtxt = 'Advection properties have not been set';
            advobj = getClassObj(mobj,'Inputs','Advection',msgtxt);
            if isempty(advobj), return; end
            msgtxt = 'Water Level properties have not been set';
            wlvobj  = getClassObj(mobj,'Inputs','WaterLevels',msgtxt);
            if isempty(wlvobj), return; end

            %set up the initial plot of tidal pumping discharge  
            sfact = 10;   %set slider range as factor of max river flow
            Qin = max([obj(:).RiverFlow]);
            if Qin==0, return; end
            Qr = [Qin/sfact,Qin*sfact];
            [qtp0,reachChIDs] = getTidalPumpingFlows(obj,mobj,1);
            if isempty(qtp0), return; end
            
            %get distances from mouth
            rchobj = getClassObj(mobj,'Inputs','Reach');
            xi = [rchobj(:).CumulativeLength];
            [xi,idxi] = sort(xi);   %sort data into distance order
            qtp0 = qtp0(idxi);
            
            %remove any existing plot and generate new plot on tab
            ht = findobj(src,'Type','axes');            
            delete(ht);
            ax = axes('Parent',src,'Tag','PlotFigAxes');
            ax.Position = [0.16,0.18,0.65,0.75]; %make space for slider bar            
            setSlideControl(obj,src,Qr(1),Qr(end),Qin);
            legtxt = sprintf('TP for Q = %.1f m^3/s',Qin);
            
            %create base plot
            stem(ax,xi,qtp0,'xb','DisplayName','Current settings')
            hold on
%               plot(ax,xi,qtp,'Color',mcolor('orange'),'DisplayName',legtxt)
              s1 = stem(ax,xi,qtp0,'Color',mcolor('orange'),...
                      'DisplayName',legtxt,...
                      'ButtonDownFcn',@(src,evt)setDataTip(obj,src,evt));
%               s1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
              s1.UserData = struct('id',[0;reachChIDs],'order',idxi);
            hold off
            xlabel('Distance from mouth (m)')
            ylabel('Discharge (m^3/s)')
            title('Tidal pumping based on current settings')
            legend
            ax.UserData = mobj;
        end
    end
    
%%
    methods (Access=private)
        function [obj,idx] = editRiver(obj,eleobj,flag)
            %edit the properties of a river instance            
            numRivers = length(obj);
            riversList = cell(numRivers,1);
            elenames = getEleProp(eleobj,'EleName');
            for i=1:numRivers
                channelid = num2str(obj(i).ChannelID);
                idele = [eleobj(:).EleID]==obj(i).ChannelID;
                riversList{i} = sprintf('%s: %s',channelid,elenames{idele});
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
        function [qtp,reachChIDs] = getTidalPumpingFlows(~,mobj,qfact)
            %allow the user to test the influence of varying river discharge
            %on tidal pumping
            asmobj.Time = 0;   %initialise water levels for t=0
            WaterLevels.setWaterLevels(mobj,asmobj);
            %initialise advection network and reach properties
            advobj = getClassObj(mobj,'Inputs','Advection');
            riverGraph = Advection.initialiseRiverGraph(mobj);
            FlowIn = advobj.RiverIn;

            %get reach based tidal pumping discharge for user Qin
            exchIn = FlowIn*qfact;
            riverGraph = rescale_graph(riverGraph,exchIn,true);
            advobj.RiverGraph = riverGraph;
            setClassObj(mobj,'Inputs','Advection',advobj);
            Reach.setReach(mobj,true); %sets IDs and properties
            
            g_landward = Reach.getReachGraph(mobj);
            reachChIDs = g_landward.Nodes.EleID(2:end);
            [tp,tpM] = Advection.getTidalPumpingDischarge(mobj);  
            qtp = [tpM(1);tp];
        end
%%
        function hm = setSlideControl(obj,hfig,qmin,qmax,qin)
            %intialise slider to set different Q values 
            invar = struct('sval',[],'smin',[],'smax',[],'size', [],...
                           'callback','','userdata',[],'position',[],...
                           'stxext','','butxt','','butcback','');            
            invar.sval = qin;        %initial value for slider 
            invar.smin = qmin;       %minimum slider value
            invar.smax = qmax;       %maximum slider value
            invar.size = [0.1,0.1];  %stepsize for slider
            invar.callback = @(src,evt)updateTPplot(obj,src,evt); %callback function for slider to use
            invar.userdata = qin;  %user data passed to widget
            invar.position = [0.15,0.005,0.45,0.04]; %position of slider
            invar.stext = 'River discharge = ';   %text to display with slider value, if included          
            hm = setfigslider(hfig,invar);   
        end 
%%
        function setDataTip(~,src,evt)
            %add data tip with EleID, distance and discharge at selected point
            dt = findobj(src,'Type','datatip');
            if isempty(dt) 
                dt = datatip(src,evt.IntersectionPoint(1),evt.IntersectionPoint(2));
                if length(src.DataTipTemplate.DataTipRows)==2
                    row = dataTipTextRow('Ele ID: ',src.UserData.id);
                    src.DataTipTemplate.DataTipRows(end+1) = row;
                end
            end
            src.DataTipTemplate.DataTipRows(1).Label = 'Upstream x: ';
            src.DataTipTemplate.DataTipRows(2).Label = 'Qtp: ';           
        end
%%
        function updateTPplot(obj,src,~)
            %use the updated slider value to adjust the TidalPumping plot
            stxt = findobj(src.Parent,'Tag','figsliderval');
            Q = round(src.Value);
            stxt.String = num2str(Q);     %update slider text
            sldui = findobj(src.Parent,'Tag','figslider');
            qfact = Q/sldui.UserData;
            
            %figure axes and update plot
            ax = findobj(src.Parent,'Type','axes');
            [qtp,~] = getTidalPumpingFlows(obj,ax.UserData,qfact);             
            hline = findobj(ax,'-regexp','DisplayName','TP for Q');
            xi = hline.XData;
            s1user = hline.UserData;
            qtp = qtp(s1user.order);
            delete(hline)
            legtxt = sprintf('TP for Q = %d m^3/s',Q);
            hold on
            s1 = stem(ax,xi,qtp,'Color',mcolor('orange'),...
                      'DisplayName',legtxt,...
                      'ButtonDownFcn',@(src,evt)setDataTip(obj,src,evt));
            s1.UserData = s1user;
            hold off
            legend
        end
%%        
        function [rownames,colnames,userdata] = riverProperties(obj,mobj)
            %define river Properties  
            advobj = getClassObj(mobj,'Inputs','Advection');
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
                if isfield(advobj,'RiverIn') && ~isempty(advobj.RiverIn)                            
                    userdata{2,k} = num2str(advobj.RiverIn(obj(j).ChannelID));
                end  
                userdata{3,k} = num2str(obj(j).RiverRho); 
            end
        end
        
    end
end