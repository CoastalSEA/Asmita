classdef Interventions < matlab.mixin.Copyable
    %class to define and update changes that are imposed on the system
    %----------------------------------------------------------------------
    % NOTES
    % When interventions added for first time, selectIntElement is called
    % and creates an instance for each element in the model. Each instance
    % has a dummy struct with the start year of the model and zeros.
    %
    % When using Intervention of tidal flats and other elements that can be
    % defined as water or sediment volumes, it is advisable to use water
    % volumes (i.e. n>0). However, the code should work for both cases.
    %
    % Volumes and surface areas that are <=0 are set to a value of 999 to
    % avoid matrices becoming illconditioned.
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
        Year           %values for each elemet are stored in a vector array
        VolumeChange 
        SurfaceAreaChange 
        isNonErodible 
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
            %select element and assign interventions using table dialogue
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
            
            %for each element sort the specified changes into a vectors of
            %unique years and volume/area changes
            [uyrs,vals,isnero] = sortInterventions(obj);
            nintele = length(obj);  %number of elements with interventions
            for i=1:nintele
                obj(i).transVolChange = 0;  %initialise transient properties
                obj(i).transAreaChange = 0;
            end
            setClassObj(mobj,'Inputs','Interventions',obj);
            
            %store of changes set up during model initialisation
            asmobj = getClassObj(mobj,'Inputs','ASM_model');
            asmobj.UniqueYears = uyrs; 
            asmobj.AnnualChange = vals;
            asmobj.isFixed = isnero;
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
            obj = blankInstance(obj,eleid,mobj);
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
            obj(eleid).Year = data{1};
            obj(eleid).VolumeChange = data{2};
            obj(eleid).SurfaceAreaChange = data{3};
            obj(eleid).isNonErodible = data{4};
            setClassObj(mobj,'Inputs','Interventions',obj);
        end     

%%
        function [idx,ok] = setAnnualChange(mobj,robj)
            %called at each time step to determine whether there any
            %changes within the time step.
            idx = [];   ok = 1; 
            obj  = getClassObj(mobj,'Inputs','Interventions');
            if isempty(obj), return; end
            eleobj  = getClassObj(mobj,'Inputs','Element');

            nele = length(eleobj);
            y2s = mobj.Constants.y2s;
            year = robj.DateTime;
            dt = robj.delta;
            asmobj = getClassObj(mobj,'Inputs','ASM_model');
            uyrs = asmobj.UniqueYears*y2s; 
            vals = asmobj.AnnualChange;
            isfix = asmobj.isFixed;
            idt = [];
            for j=1:length(uyrs)
                if year<=uyrs(j) && year+dt>uyrs(j) %ensures only assigned 
                    idt = j;                        %within one time interval
                end
            end
            %if no change in given year, no update needed
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            if isempty(idt) || ~rncobj.IncInterventions
                for i=1:nele
                    obj(i).transVolChange = 0;
                    obj(i).transAreaChange = 0;
                end
                return;
            end
            dV = vals(:,idt,1);
            dS =  vals(:,idt,2);
            idx = dV>0 | dS>0;
            n = getEleProp(eleobj,'TransportCoeff');
            sgn = sign(n);
            %update transient volume and area properties
            %convention is that +ve is an increase in water volume
            %so subtract for elements defined as sediment volumes
            for i=1:nele
                ok = checkInterventions(obj,eleobj,sgn(i)*dV(i),dS(i),i);
                if ok<1, return; end
                eleobj(i).MovingVolume = setVarChange(obj,eleobj(i).MovingVolume,sgn(i)*dV(i));
                eleobj(i).FixedVolume = eleobj(i).FixedVolume+sgn(i)*dV(i); %NB fixed volume can be -ve

                if dS(i)==999
                    eleobj(i).SurfaceArea = dS(i);
                else
                    eleobj(i).SurfaceArea = setVarChange(obj,eleobj(i).SurfaceArea,dS(i));
                    eleobj(i).EqSurfaceArea = setVarChange(obj,eleobj(i).EqSurfaceArea,dS(i));
                    obj(i).transVolChange = dV(i);
                    obj(i).transAreaChange = dS(i);
                end

                if isfix(i,idt)
                    %store the cumulative changes that are non-erodible to
                    %adjust EqVolume each time it is reset based on prism
                    %tidal range or drift rate. This depends on what is
                    %changing, equilibrium condition and type of element.
                    [isok,~] = checkFixedConditions(obj,mobj,dV(i),dS(i),i);
                    eleobj(i).eqFixedInts(1) = eleobj(i).eqFixedInts(1)+sgn(i)*dV(i)*isok;
                    %eleobj(i).eqFixedInts(2) = eleobj(i).eqFixedInts(2)+dS(i)*isok;
                end               
            end
            %assign updated instances of Element and Interventions
            setClassObj(mobj,'Inputs','Interventions',obj);
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
            rnpobj = getClassObj(mobj,'Inputs','RunProperties');
            
            ht = findobj(src,'Type','Axes');
            delete(ht);
            
            answer = questdlg('Select to plot Element, All or Depth','Interventions',...
                                            'Element','All','Depth','Element');
            if strcmp(answer,'Element')
                %check which elements have interventions and allow user selection
                [idx,legtxt,ok] = selectInterventionSet(intobj,mobj);
                if ok<1, return; end
                tim = intobj(idx).Year;
                Vol = cumsum(intobj(idx).VolumeChange);
                Surf = cumsum(intobj(idx).SurfaceAreaChange);
                isfix = intobj(idx).isNonErodible;
            elseif  strcmp(answer,'All')
                [tim,vals,isnero] = sortInterventions(intobj);
                Vol = sum(vals(:,:,1),1)';
                Surf = sum(vals(:,:,2),1)';
                isfix = any(isnero,1)';
                legtxt  = 'Cumulative change for All elements';
            else
                eleobj  = getClassObj(mobj,'Inputs','Element');
                if isempty(eleobj), return; end
                V0 = getEleProp(eleobj,'InitialVolume');
                S0 = getEleProp(eleobj,'InitialSurfaceArea');
                [tim,vals,~] = sortInterventions(intobj);
                
                Vol = V0+cumsum(vals(:,:,1),2);
                Surf = S0+cumsum(vals(:,:,2),2);
                
                if tim(1)>rnpobj.StartYear
                    tim = [rnpobj.StartYear;tim];
                    Vol = [V0,Vol];
                    Surf = [S0,Surf];
                end
                Depth = Vol./Surf;
                legtxt  = getEleProp(eleobj,'EleName');
                depthPlot(intobj,src,tim,Depth,legtxt);
                return
            end

            %if the run properties have been defined use these to adjust
            %time range of plot            
            if ~isempty(rnpobj) && length(Vol)>5  %5 used as switch to stair plot in tabPlot
                t0 = rnpobj.StartYear;
                tN = t0+rnpobj.TimeStep*rnpobj.NumSteps;
                if t0<tim(1)
                    tim = [t0;tim(1);tim];
                    Vol = [0;0;Vol];
                    Surf = [0;0;Surf];
                    isfix = [0;0;isfix];
                end
                %
                if tN>tim(end)
                    tim = [tim;tN];
                    Vol = [Vol;Vol(end)];
                    Surf = [Surf;Surf(end)];
                    isfix = [isfix;isfix(end)];
                end
            end
            tabPlot(intobj,src,tim,Vol,Surf,isfix,legtxt);
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
            obj(idx).VolumeChange = -1*obj(idx).VolumeChange;
            obj(idx).SurfaceAreaChange = -1*obj(idx).SurfaceAreaChange;
        end
    end
%%
    methods (Static,Hidden)
        function setNewIntervention(mobj)
            %initialise an empty instance of Interevention (used in asm_oo2mui)
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
            obj = blankInstance(obj,nint,mobj);
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
            userdata = {obj(eleid).Year,...
                        obj(eleid).VolumeChange,...
                        obj(eleid).SurfaceAreaChange,...
                        obj(eleid).isNonErodible};
            colnames = {'Year','Volume','SurfaceArea','NonErodible'};
            nyr = length(obj(eleid).Year);
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
            
            obj(eleid).Year = newtable.Year;
            obj(eleid).VolumeChange = newtable.Volume;
            obj(eleid).SurfaceAreaChange = newtable.SurfaceArea;
            obj(eleid).isNonErodible = newtable.NonErodible;
        end

%%
        function [idx,legtxt,ok] = selectInterventionSet(obj,mobj)
            %determine which elements have interventions defined and allow
            %user to select from subset
            idx = []; legtxt = ''; ok = 0;
            eleobj = getClassObj(mobj,'Inputs','Element');
            styear = startYear(obj,mobj);
            count = 0;
            intid = [];
            for i=1:length(obj)                
                if any(obj(i).VolumeChange~=0) || ...
                          any(obj(i).SurfaceAreaChange~=0) || ...
                                  any(obj(i).Year>styear)
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
        function tabPlot(~,src,tim,Vol,Surf,isfix,legtxt)
            %plot cumulative changes based on user specified changes at given t
            offset = 0.1;  %offset from max value to provide margin at top
            vmax = max(Vol); smax = max(Surf); 
            fVol = Vol.*logical(isfix); fVol(fVol==0) = NaN;
            fSurf = Surf.*logical(isfix); fSurf(fSurf==0) = NaN;
            axes('Parent',src); 
            
            s1 = subplot(2,1,1); 
            if length(Vol)<5
                bar(tim,Vol,'DisplayName',legtxt);
                hold on 
                bar(tim,fVol);
                hold off
            else
                stairs(tim,Vol,'DisplayName',legtxt,'LineWidth',1);
                hold on
                stairs(tim,fVol,'LineWidth',1);
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
            legend('Erodible','Non-erodible','Location','northwest')
            %
            s2 =  subplot(2,1,2);
            if length(Surf)<5
                bar(tim,Surf,'DisplayName',legtxt);
                hold on 
                bar(tim,fSurf);
                hold off
            else
                stairs(tim,Surf,'DisplayName',legtxt,'LineWidth',1);
                hold on
                stairs(tim,fSurf,'LineWidth',1);
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

%%
        function depthPlot(~,src,tim,Depth,legtxt)
            %plot hydraulic depths for all elements including interventions
            ax = axes('Parent',src);
            dmin = min(Depth,[],'All');
            plot(ax,tim,Depth(1,:))
            hold on
            for i=2:size(Depth,1)
                plot(ax,tim,Depth(i,:))
            end
            if dmin<0
                hp = plot(xlim,[0,0],'-.r');
                hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            end
            hold off
            xlabel('Year');
            ylabel('Depth (m)');
            legend(legtxt,'Location','best')
        end

%%
        function [uyrs,vals,isnero] = sortInterventions(obj)
            %for each element sort the specified changes into vectors of
            %unique years and volume/area changes
            nintele = length(obj);  %number of elements with interventions
            uyrs = [];
            for i=1:nintele
                uyrs = cat(1,uyrs,obj(i).Year);
            end
            uyrs = unique(uyrs);  %returns sorted order
            nyr  = length(uyrs);
            vals = zeros(nintele,nyr,2);
            isnero = zeros(nintele,nyr);
            for i=1:nintele
                year = obj(i).Year;
                vols = obj(i).VolumeChange;
                area = obj(i).SurfaceAreaChange;
                flag = obj(i).isNonErodible;
                for j=1:nyr
                    if any(year==uyrs(j))
                        vals(i,j,1) = sum(vols(year==uyrs(j)));
                        vals(i,j,2) = sum(area(year==uyrs(j)));
                        isnero(i,j) = flag(year==uyrs(j));
                    end
                end
            end
        end

%%
        function ok = checkInterventions(~,eleobj,dV,dS,idx)
            %check that the intervention is not bigger than the element
            errtxt = []; ok = 1;
            if eleobj(idx).MovingVolume+dV<0
                errtxt = sprintf('Intervention in Element %d results in negative volume',idx);
            elseif eleobj(idx).SurfaceArea+dS<0
                errtxt = sprintf('Intervention in Element %d results in negative surface area (Sa)',idx);
            elseif eleobj(idx).EqSurfaceArea+dS<0
                errtxt = sprintf('Intervention in Element %d results in negative surface area (Se)',idx);
            end
        
            if ~isempty(errtxt)
                warndlg(sprintf('%s\nRun aborted in Interventions.setAnnualChange',errtxt))
                ok = 0; return;
            end
        end

%%
        function Var = setVarChange(~,Vari,dVar)
            %set the change in volume or surface area
            if Vari==999 && dVar>0        %reset if change restores element
                Vari = 0;                 %remove 999 offset
            end
            
            Var = Vari+dVar;              %add intervention
            if Var<999, Var = 999; end    %adjust negative or small values
        end

%%
        function [isupdate,fcase] = checkFixedConditions(~,mobj,dV,dS,i)
            %check equilibirum conditions being used to determine how fixed
            %interventions should be applied            
            eleobj  = getClassObj(mobj,'Inputs','Element');
            eletype = eleobj(i).EleType;
            ecpobj = getClassObj(mobj,'Inputs','EqCoeffParams');
            isTReq = ~logical(ecpobj.eqtype.(eletype)); %switch to true if tidal range equilibrium
            isdV = dV~=0;   
            isdS = dS~=0;
            isdVdS = isdV && isdS;
            rchtypes = mobj.GeoType(mobj.REtypes); %reach element types
            flatypes =rchtypes(~strcmp(rchtypes,'FloodDelta')); %exclude FloodDelta
            isReachType = any(strcmp(eletype,flatypes));
            isPrismType = any(strcmp(eletype,mobj.GeoType(mobj.RPtypes)));
            if (isdVdS || isdS) && isReachType && isTReq
                %dV and dS, or dS, change for elements that can belong 
                %to a reach using tidal range to define equilibrium
                isupdate = false; fcase = 1;
            elseif isdV && isReachType && isTReq
                %dV changes for elements that can belong to a reach
                %using tidal range to define equilibrium
                isupdate = true; fcase = 2;
            elseif isdV && isReachType && ~isTReq 
                %dV changes for elements that can belong to a reach
                %using tidal prism to define equilibrium
                isupdate = false; fcase = 3;
            elseif (isdVdS || isdS) && isPrismType && ~isTReq   
                %dV changes for elements with prism equilibrium that do not
                %alter prism (channel, flood delta, ebb delta) Note that
                %isdS case only applies to channel and flood delta but as
                %false this case can be included here.
                isupdate = false; fcase = 4;
            elseif isdV && isPrismType && ~isTReq
                %dV changes for elements with prism equilibrium that do not
                %alter prism (channel, flood delta, ebb delta)
                isupdate = true; fcase=5;
            else 
                isupdate = false; fcase = [];
            end
        end
        %%
        function obj = blankInstance(obj,eleid,mobj)
            %assign dummy values to an interventions struct
            obj(eleid).Year = startYear(obj,mobj);
            obj(eleid).VolumeChange = 0;
            obj(eleid).SurfaceAreaChange = 0;
            obj(eleid).isNonErodible = 0;
        end

%%
        function styear = startYear(~,mobj)
            %get start year for model if set
            rnpobj  = getClassObj(mobj,'Inputs','RunProperties');
            if isempty(rnpobj) 
                styear = 0; 
            else
                styear = rnpobj.StartYear;
            end
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

            if ~isempty(obj) && length(obj)~=length(eleobj)
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
                    obj = blankInstance(obj,i,mobj);
                end
            end
            %
            [eleid,elename,ok] = selectElement(eleobj);
        end
    end
end



