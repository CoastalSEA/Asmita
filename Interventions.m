classdef Interventions < matlab.mixin.Copyable
    %class to define and update changes that are imposed on the system
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
        ElementID      %should be a copy of Element.EleID (cross-check)
        Year = {0}     %values for each elemet are stored in a cell array
        VolumeChange = {0}
        SurfaceAreaChange = {0}
    end
    
    properties (Transient)
        transVolChange = 0
        transAreaChange = 0
    end

    methods (Access=private)
        function obj = Interventions
        end
    end
%%    
    methods (Static)
        function setInterventions(mobj)
            %select element and assign interventioons using table dialogue
            [obj,eleid,elename,ok] = Interventions.selectIntElement(mobj);
            if ok<1, return; end
            obj = intTable(obj,eleid,elename);
            setClassObj(mobj,'Inputs','Interventions',obj);
        end       
%%         
        function initialiseInterventions(mobj)
            %initialise the transient intervention properties and set the
            %years and annual change values in ASM_model
            obj  = getClassObj(mobj,'Inputs','Interventions');
            if isempty(obj), return; end
            
            %for each element sort the specified changes into an vectors of
            %unique years and volume/area changes
            nintele = length(obj);  %number of elements with interventions
            uyrs = [];
            for i=1:nintele
                uyrs = cat(1,uyrs,obj(i).Year{1});
            end     
            uyrs = unique(uyrs);  %returns sorted order
            nyr  = length(uyrs);
            vals = zeros(nintele,nyr,2);
            for i=1:nintele
                obj(i).transVolChange = 0;  %initialise transient properties
                obj(i).transAreaChange = 0;
                year = obj(i).Year{1};
                vols = obj(i).VolumeChange{1};
                area = obj(i).SurfaceAreaChange{1};
                for j=1:nyr
                    vals(i,j,1) = sum(vols(year==uyrs(j)));
                    vals(i,j,2) = sum(area(year==uyrs(j)));
                end
            end
            setClassObj(mobj,'Inputs','Interventions',obj);
            
            %store of changes set up during model initialisation
            asmobj  = getClassObj(mobj,'Inputs','ASM_model');
            asmobj.UniqueYears = uyrs; 
            asmobj.AnnualChange = vals;
            setClassObj(mobj,'Inputs','ASM_model',asmobj);
        end      
%%
        function clearEleInt(mobj,eleid)
            %clear and reset an intervention from an element
            % called from main UI
            msgtxt = 'No interventions have been defined';
            obj  = getClassObj(mobj,'Inputs','Interventions',msgtxt);
            if isempty(obj), return; end
            %
            if nargin<2
                [eleid,~,ok] = selectInterventionSet(obj,mobj);
                if ok<1, return; end
            end
            %
            obj(eleid).Year = {0};
            obj(eleid).VolumeChange = {0};
            obj(eleid).SurfaceAreaChange = {0};
            setClassObj(mobj,'Inputs','Interventions',obj);
        end
%%
        function loadInterventions(mobj)
            %prompt user for file name and setup intervention timeseries
            %each element can have an independent timeseries that is
            %loaded from a file (one element at a time)
            [obj,eleid,~,ok] = Interventions.selectIntElement(mobj);
            if ok<1, return; end
            %get the name of the file and read it
            data = readinputdata();        %see asmita_functions (different to readinputfile)
            if isempty(data), return; end  %user aborted or file not read
            
            %assign interventions to object
            obj(eleid).Year = data(1);
            obj(eleid).VolumeChange = data(2);
            obj(eleid).SurfaceAreaChange = data(3);
            setClassObj(mobj,'Inputs','Interventions',obj);
        end         
%%
        function setAnnualChange(mobj,robj)
            %called at each time step to determine whether there any
            %changes within the time step.
            intobj  = getClassObj(mobj,'Inputs','Interventions');
            if isempty(intobj), return; end
            eleobj  = getClassObj(mobj,'Inputs','Element');

            nele = length(eleobj);
            y2s = mobj.Constants.y2s;
            year = robj.DateTime;
            dt = robj.delta;
            asmobj = getClassObj(mobj,'Inputs','ASM_model');
            uyrs = asmobj.UniqueYears*y2s; 
            vals = asmobj.AnnualChange;
            idx = [];
            for j=1:length(uyrs)
                if year<=uyrs(j) && year+dt>uyrs(j) %ensures only assigned 
                    idx = j;                        %within one time interval
                end
            end
            %if no change in given year, no update needed
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            if isempty(idx) || ~rncobj.IncInterventions
                for i=1:nele
                    intobj(i).transVolChange = 0;
                    intobj(i).transAreaChange = 0;
                end
                return;
            end
            dV = vals(:,idx,1);
            dS =  vals(:,idx,2);
            n = sign(getEleProp(eleobj,'TransportCoeff'));
            %update transient volume and area properties
            %convention is that +ve is an increase in water volume
            %so subtract for elements defined as sediment volumes
            for i=1:nele
                eleobj(i).MovingVolume = eleobj(i).MovingVolume+n(i)*dV(i);
                eleobj(i).MovingSurfaceArea = eleobj(i).MovingSurfaceArea+dS(i);
                eleobj(i).FixedVolume = eleobj(i).FixedVolume+n(i)*dV(i);
                eleobj(i).FixedSurfaceArea = eleobj(i).FixedSurfaceArea+dS(i);
                intobj(i).transVolChange = dV(i);
                intobj(i).transAreaChange = dS(i);
            end
            %assign updated instances of Element and Interventions
            setClassObj(mobj,'Inputs','Interventions',intobj);
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
%%
        function prop = getIntProp(mobj,varname)
            %get intervention property and return as an element array
            prop = 0;
            obj  = getClassObj(mobj,'Inputs','Interventions');
            if isempty(obj), return; end           
            prop = [obj.(varname)]';
        end      
%%    
        function IntTabPlot(mobj,src,~)
            %plot the interventions that have been defined on the tab
            msgtxt = 'No interventions have been defined';
            intobj  = getClassObj(mobj,'Inputs','Interventions');
            if isempty(intobj), getdialog(msgtxt); return; end    
            
            ht = findobj(src,'Type','Axes');
            delete(ht);
            
            %check which elments have interventions and allow user selection
            [idx,legtxt,ok] = selectInterventionSet(intobj,mobj);
            if ok<1, return; end

            tim = intobj(idx).Year{:};
            Vol = cumsum(intobj(idx).VolumeChange{:});
            Surf = cumsum(intobj(idx).SurfaceAreaChange{:});
            
            %if the run properties have been defined use these to adjust
            %time range of plot
            rnpobj = getClassObj(mobj,'Inputs','RunProperties');
            if ~isempty(rnpobj) && length(Vol)>5  %5 used as switch to stair plot in tabPlot
                t0 = rnpobj.StartYear;
                tN = t0+rnpobj.TimeStep*rnpobj.NumSteps;
                if t0<tim(1)
                    tim = [t0;tim(1);tim];
                    Vol = [0;0;Vol];
                    Surf = [0;0;Surf];
                end
                %
                if tN>tim(end)
                    tim = [tim;tN];
                    Vol = [Vol;Vol(end)];
                    Surf = [Surf;Surf(end)];
                end
            end
            tabPlot(intobj,src,tim,Vol,Surf,legtxt);
        end       
%%
        function setIntChSign(mobj)
            %change the sign of interventions for selected element
            %Useful if wrong convention has been used on data entry.
            % Convention: Loss of water volume and/or plan area
            % should be negative
            msgtxt = 'No interventions have been defined';
            obj  = getClassObj(mobj,'Inputs','Interventions',msgtxt);
            if isempty(obj), return; end 
            %check which elments have interventions and allow user selection
            [idx,~,ok] = selectInterventionSet(obj,mobj);
            if ok<1, return; end
            %change sign of interventions for selected element
            obj(idx).VolumeChange = {-1*obj(idx).VolumeChange{1}};
            obj(idx).SurfaceAreaChange = {-1*obj(idx).SurfaceAreaChange{1}};
        end
    end
%%
    methods (Static,Hidden)
        function setNewIntervention(mobj)
            %initialise an empty instance ofInterevention (used in asm_oo2mui)
            obj = Interventions;
            setClassObj(mobj,'Inputs','Interventions',obj);
        end
    end
%%
    methods
        function addIntEle(obj,mobj,elobj)
            %add an intervention to the array of intervention objects when
            %an element is added - called from Element.addElement
            nint = length(obj)+1;
            if nint~=length(elobj)
                warndlg('Error in Interventions.addIntEle')
                return
            end
            obj(nint) = Interventions;
            obj(nint).ElementID = elobj(nint).EleID;
            setClassObj(mobj,'Inputs','Interventions',obj);
        end        
%%
        function delIntEle(obj,mobj,eleid)
            %delete an intervention from an element
            % used in Element.delElement
            obj(eleid) = [];
            setClassObj(mobj,'Inputs','Interventions',obj);
        end        
    end
%%
    methods (Access=private)
        function obj = intTable(obj,eleid,elename) 
            %generate UI table for user to edit and add to element interventions
            userdata = {obj(eleid).Year{1},...
                        obj(eleid).VolumeChange{1},...
                        obj(eleid).SurfaceAreaChange{1}};
            colnames = {'Year','Volume','SurfaceArea'};
            nyr = length(obj(eleid).Year{1});
            
            %create row index
            rownames{nyr,1} = '';
            for i=1:nyr, rownames{i}=num2str(i); end

            title = 'Interventions';
            header = sprintf('Enter interventions for %s',elename);
            oldtable = table(userdata{:},'VariableNames', colnames); 
            but.Text = {'Save','Add','Cancel'}; %labels for tab button definition
            newtable = tablefigureUI(title,header,oldtable,true,but,[0.1,0.6]);
            if isempty(newtable), return; end 
            newtable = sortrows(newtable);
            %remove row if year has not been defined
            newtable = rmmissing(newtable,'DataVariables',{'Year'});
            %if user has not edited all NaN values amend to zeros
            newtable = fillmissing(newtable,'constant',0);
            
            obj(eleid).Year = {newtable.Year};
            obj(eleid).VolumeChange = {newtable.Volume};
            obj(eleid).SurfaceAreaChange = {newtable.SurfaceArea};
        end
%%
        function [idx,legtxt,ok] = selectInterventionSet(obj,mobj)
            %determine which elements have interventions defined and allow
            %user to select from subset
            idx = []; legtxt = ''; ok = 0;
            eleobj = getClassObj(mobj,'Inputs','Element');
            count = 0;
            intid = [];
            for i=1:length(obj)                
                if any(obj(i).VolumeChange{1}~=0) || ...
                                        any(obj(i).SurfaceAreaChange{1}~=0) 
                    count = count+1;
                    intid(count) = i;                     %#ok<AGROW>
                    eleid = obj(i).ElementID;
                    elename{count,1} = eleobj(eleid).EleName;  %#ok<AGROW>
                end               
            end
            %select element to use
            if isempty(intid)
                    warndlg('No intervention values defined')                    
                    return;
            elseif length(intid)>1
                [select, ok] = listdlg('Name','Interventions', ...
                    'PromptString','Select Element with Interventions', ...
                    'SelectionMode','single', ...
                    'ListString',elename);   
                if ok==0, return; end
            else
                select = 1; ok = 1;
            end 
            idx = intid(select);
            eleid = obj(idx).ElementID;
            legtxt = sprintf('Cumulative change for element %d: %s',eleid,char(elename{select}));
        end
%%
        function tabPlot(~,src,tim,Vol,Surf,legtxt)
            %plot cumulative changes based on user specified changes at given t
            offset = 0.1;  %offset from max value to provide margin at top
            vmax = max(Vol); smax = max(Surf);
            axes('Parent',src);            
            s1 = subplot(2,1,1); 
            if length(Vol)<5
                bar(tim,Vol,'DisplayName',legtxt);
            else
                stairs(tim,Vol,'DisplayName',legtxt,'LineWidth',1);
                hold on
                p1 = plot([tim(1),tim(end)],[Vol(1),Vol(end)],'.','MarkerSize',10);
                p1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
                hold off
            end
            title(legtxt);
            ylabel('Volume change (m^3)');
            if vmax<0 && vmax<s1.YLim(2)
                s1.YLim(1) = vmax+offset*vmax;
            else
                s1.YLim(2) = vmax+offset*vmax;
            end
            s1.XLim(2) = tim(end)+1;
            %
            s2 =  subplot(2,1,2);
            if length(Surf)<5
                bar(tim,Surf,'DisplayName',legtxt);
            else
                stairs(tim,Surf,'DisplayName',legtxt,'LineWidth',1);
                hold on
                p2 = plot([tim(1),tim(end)],[Surf(1),Surf(end)],'.','MarkerSize',10);
                p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
                hold off
            end
            xlabel('Year');
            ylabel('Area change (m^2)');
            if smax<=0 && smax<s2.YLim(2)
                s2.YLim(1) = smax+offset*smax;
            else
                s2.YLim(2) = smax+offset*smax;
            end
            s2.XLim(2) = tim(end)+1;
        end
    end
    
%%    
    methods (Access=private, Static)      
        function [obj,eleid,elename,ok] = selectIntElement(mobj)
            %select an element to use for interventions input/modification
            % creates array of class instances if not already instantiated
            eleid = []; elename = []; ok = 0;
            msgtxt = 'No Elements have been defined';            
            eleobj  = getClassObj(mobj,'Inputs','Element',msgtxt);
            if isempty(eleobj), obj = eleobj; return; end

            obj  = getClassObj(mobj,'Inputs','Interventions');

            if length(obj)~=length(eleobj)
                %interventions dims do not match element dims
                msgerr = sprintf('Element dimensions do not match\nInterventions have been cleared\nYou may need to input intervention data again.');
                warndlg(msgerr)
                obj = [];
            end

            if isempty(obj)    %no inverventions defined
                clear obj
                nele = length(eleobj);
                obj(nele,1) = Interventions;
                for i=1:nele
                    obj(i).ElementID = eleobj(i).EleID;
                end
            end
            %
            [eleid,elename,ok] = selectElement(eleobj);
        end
    end
end



