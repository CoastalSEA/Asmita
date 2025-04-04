classdef ASM_Plots < muiPlots
%
%-------class help---------------------------------------------------------
% NAME
%   ASM_Plots.m
% PURPOSE
%   Class to implement the generation of a range of plot types in
%   CoastalTools
% NOTES
%   inherits from muiPlots
% SEE ALSO
%   called from ASM_PlotsUI.m, which defines the selection and settings in
%   properties UIselection and UIsettings. Uses muiCatalogue to access
%   data. 
%
% Author: Ian Townend
% CoastalSEA (c)July 2021
%--------------------------------------------------------------------------
%
    properties (Transient)
        %inherits the following properties from muiPlots
        %Plot            %struct array for:
                         %   FigNum - index to figures created
                         %   CurrentFig - handle to current figure
                         %   Order - struct that defines variable order for
                         %           plot type options (selection held in Order)

        %ModelMovie      %store animation in case user wants to save
        %UIsel           %structure for the variable selection made in the UI
        %UIset           %structure for the plot settings made in the UI
        %Data            %data to use in plot (x,y,z)
        %TickLabels      %struct for XYZ tick labels
        %AxisLabels      %struct for XYZ axis labels
        %Legend          %Legend text
        %MetaData        %text summary of primary variable selection
        %Title           %Title text
        %Order           %order of variables for selected plot type
        %idxfig          %figure number of the current figure
    end
%%
    methods 
        function obj = ASM_Plots
            %types of plot avaiable based on number of dimensions
            obj.Plot.FigNum = [];
            obj.Plot.CurrentFig = [];  
            obj.Plot.Order = ASM_Plots.setVarOrder;
        end      
    end
%%    
    methods (Static)
        function getPlot(gobj,mobj)
            %get existing instance or create new class instance
            if isa(mobj.mUI.Plots,'ASM_Plots') 
                obj = mobj.mUI.Plots;    %get existing instance          
                clearPreviousPlotData(obj);
            else
                obj = ASM_Plots;                   %create new instance
            end

            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;

            %set the variable order for selected plot type
            obj.Order = obj.Plot.Order.(obj.UIset.callTab);            
            
            %modify the legend text to include case and variable only
            legformat.idx = [0 1 0];  %format as Case(Element)Variable
            legformat.text = {};
            
            dtype = 'array';
            switch obj.UIset.callTab
                case 'Time'
                    obj = addSelection(obj,mobj);                   
                    legformat.text = {obj.UIsel(1).Element};
                case {'Distance','Network'}
                    if strcmp(obj.UIset.callTab,'Distance') && isSnapShot(obj,2)
                        %get the current slider value if a distance snap shot plot
                        obj.UIset.RunTime = gobj.TabContent(2).Selections{4}.Value;
                        %format as Case(Time)Variable
                        legformat.text = {num2str(obj.UIset.RunTime,'%.1f')};
                    else
                        legformat.idx = [1 0 1]; %format as Case(Variable)
                    end
                    obj = addSelection(obj,mobj);                        
                case {'2D','2DT'}
                    legformat = repmat(legformat,1,2);
                    legformat(1).text = obj.UIsel(1).dims(2).value;
                    legformat(2).text = obj.UIsel(2).dims(2).value;
                case {'3D','3DT'}                    
                    legformat = repmat(legformat,1,3);
                    for i=1:3
                        legformat(i).text = obj.UIsel(i).dims(2).value;
                    end
            end
            
            %get the data to be used in the plot 
            if isfield(obj.UIset,'AllEle') && obj.UIset.AllEle
                multiElePlot(obj,mobj);
                return;
            else
                ok = getPlotData(obj,mobj.Cases,dtype,legformat);
            end
            if ok<1, return; end %data not found
            
            if strcmp({'Network'},obj.UIset.callTab)
                %add the digraph data for the selected case
                obj = getNetworkData(obj,mobj);
            elseif strcmp({'Distance'},obj.UIset.callTab)
                obj = setDistancePlotType(obj,mobj);
            else
                % if contains(obj.UIset.callTab,{'2D','2DT','3D','3DT'})
                %   obj = setAxisTicks(obj);                               %no longer needed?
                % end
                %distance and network plots needs all elements so dimensions do not match
                isvalid = checkdimensions(obj);
                if ~isvalid, return; end
            end

            if ~isempty(obj.UIset.Type) && strcmp(obj.UIset.Type.String,'User')
                user_plot(obj,mobj);  %pass control to user function
            else
                %generate the plot
                setPlot(obj,mobj);  %uses muiPlots method
            end
        end
     end
%%   
    methods (Access=protected)       
        function callPlotType(obj)
            %call the function specific to the selected plot type           
            switch obj.UIset.callTab        %call function based on Tab
                case {'Time','2D'}                 
                    switch obj.UIset.callButton %and Tab Button used
                        case 'New'              %create new 2D plot
                            new2Dplot(obj);
                            if strcmp(obj.UIset.callTab,'Time')
                                addPivotButton(obj);
                            end
                        case 'Add'              %add variable to 2D plot
                            if strcmp(obj.UIset.Type,'bar')
                                addBarplot(obj);
                            else
                                add2Dplot(obj);
                            end
                        case 'Delete'           %delete variable from 2D plot
                            del2Dplot(obj);
                    end
                case 'Distance'
                    setDistancePlot(obj)
                case 'Network'
                    setNetworkPlot(obj)
                case '3D'
                    new3Dplot(obj);
                    adjustAxisTicks(obj,gca);  %adjust tick labels  if defined
                case {'2DT','3DT'}
                    newAnimation(obj);
                otherwise
                    warndlg('Could not find plot option in callPlotType');
            end
        end
%%
        function obj = addSelection(obj,mobj)
            %add a selection that is included by default (eg time or
            %direction)
            obj.UIset.Polar = false; %checked in muiPlots.newAnimation
            obj.UIsel(2) = obj.UIsel(1);
            obj.UIsel(1).xyz = [true false];
            obj.UIsel(2).xyz = ~obj.UIsel(1).xyz; %invert logical array
            
            %define time as an input
            obj.UIsel(2).variable = 1;
            obj.UIsel(2).property = 2;
            
            %assign time if subselection required using RowNames assigned
            %to obj.UIset.RunTime in getPlot above
            trange = obj.UIsel(1).dims(1).value;            
            obj.UIsel(2).range = trange;
            obj.UIsel(2).desc = 'Time';  
            obj.UIsel(2).dims = struct('name','','value',[]); 
            if strcmp(obj.UIset.callTab,'Distance') && isSnapShot(obj,2)
                runtime = datetime(0,1,1)+years(obj.UIset.RunTime);
                obj.UIsel(1).dims(1).value = runtime;  
            end
            
            %assign Element selection to dimension property
            if strcmp(obj.UIset.callTab,'Distance') || ...
                                      strcmp(obj.UIset.callTab,'Network')
                caseDef = getRunParams(obj,mobj,1);
                eleobj = caseDef.Element;
                eleprop = getEleProp(eleobj,'EleName');  %use Element table
                if obj.UIsel(1).dataset==2               %use Reach table  
                    idr = unique(getEleProp(eleobj,'ReachID'));
                    eleprop = eleprop(idr(idr>0));
                end
                obj.UIsel(1).dims(2).name = 'EleName';
                %add selection of all elements to Y plot variable
                obj.UIsel(1).dims(2).value =  eleprop;
                obj.UIsel(1).xyz = [true false false];
                obj.UIsel(2).xyz = [false true false];
                %add additional dimension for the elements
                obj.UIsel(3) = obj.UIsel(1);                
                obj.UIsel(3).range = var2range([eleprop(1),eleprop(end)]);
                obj.UIsel(3).property = 3;
                obj.UIsel(2).xyz = [false false true];
            else
                obj.UIsel(1).dims(2).name = 'EleName';
                obj.UIsel(1).dims(2).value = {obj.UIsel(1).Element};
                seltxt = obj.UIsel(1).desc;
                obj.UIsel(1).desc = replace(seltxt,'Element',obj.UIsel(1).dims(2).value{1});
            end
        end
%%
%  updates to handle lists in muiPlots means that this is not needed 
%          function obj = setAxisTicks(obj)
%             %check if any axis selection uses a character variable
%             if iscell(obj.Data.X)             
%                 obj.TickLabels.XTick = 1:length(obj.Data.X);
%                 obj.TickLabels.XTickLabel = obj.Data.X;
%                 obj.Data.X =1:length(obj.Data.X);
%             end
%             %
%             if isfield(obj.Data,'Y') && iscell(obj.Data.Y)    
%                 obj.TickLabels.YTick = 1:length(obj.Data.Y);
%                 obj.TickLabels.YTickLabel = obj.Data.Y;
%                 obj.Data.Y =1:length(obj.Data.Y);
%             end
%         end
%%
        function addPivotButton(obj)
            %add button to allow user to set pivot year
            hf = obj.Plot.CurrentFig;
            pvtcb =  @(src,evdat)asmita_plot_pivot_year(src,evdat);
            uicontrol('Parent',hf,'Style','pushbutton',...
                'String','pivot','Tag','FigButton',...
                'TooltipString','Set pivot year',...
                'Units','normalized','Position',[0.88 0.95 0.10 0.044],...
                'Callback',pvtcb);           
        end
%%
%--------------------------------------------------------------------------
% Code to handle distance plots
%--------------------------------------------------------------------------
        function setDistancePlot(obj)
            %generate new distance plot in figure  
            obj.UIset.callTab = '2DT'; %use default tab name in muiPlots
            %obj = setAxisTicks(obj); %no longer needed?
            %generate an animation of plot type is a line/point type
            if ~isSnapShot(obj,1)
                %convert from animation to 3D assignment
                obj.Data.('Z') = obj.Data.Y; 
                obj.Data.Y = obj.Data.T;
                obj.AxisLabels.Y = obj.AxisLabels.T;
                obj.Data = rmfield(obj.Data,{'T'});
                new3Dplot(obj);                
            elseif obj.UIset.Animate
                %plot selected lines/points type as an animation 
                if ~isdatetime(obj.Data.T)
                    obj.Data.T = round(obj.Data.T*10)/10;
                end
                newAnimation(obj);
            else
                %plot single snap shot
                switch obj.UIset.callButton %and Tab Button used
                    case 'Add'              %add variable to 2D plot
                        if strcmp(obj.UIset.Type,'bar')
                            addBarplot(obj);
                        else
                            add2Dplot(obj);
                        end
                    case 'Delete'           %delete variable from 2D plot
                        del2Dplot(obj);
                    otherwise               %create new 2D plot
                        new2Dplot(obj);
                end
            end
            adjustAxisTicks(obj,gca);  %adjust tick labels  if defined
        end
%%
        function obj = setDistancePlotType(obj,mobj)
            %select from x axis using elements or distances
            select = questdlg('Select plot type:','Distance plot',...
                                    'Distance','Element','Element');
            if strcmp(select,'Distance')
                obj = getDistanceData(obj,mobj);
            end
        end
%%
        function obj = getDistanceData(obj,mobj)
            %get the distance for each reach and associated variable values
            caseDef = getRunParams(obj,mobj,1);
            ismidpt = true; %returns lengths at mid-point of each element
            [reachlength,g] = Reach.getReachLengths(caseDef,ismidpt);
            yele = obj.Data.Y;
            if size(yele,2)~=length(reachlength)
                %change variable from element to reach values
                obj.Data.Y = Reach.getReachEleVar(g,yele);
            end
            %
            %if length(obj.Data.X)~=length(reachlength) %not sure what this
            %obj.Data.X  = reachlength;                 %is meant to trap????
            %end   
            obj.Data.X  = reachlength;
            obj.Title = replace(obj.Title,'Element Name','Distance');
            obj.AxisLabels.X = 'Distance (m)';
        end
%%
        function setNetworkPlot(obj)
            %call the animation method for digraphs in muiPlots
            obj.UIset.callTab = 'Digraph'; %use default tab name in muiPlots
            if ~isdatetime(obj.Data.T)
                obj.Data.T = round(obj.Data.T*10)/10;
            end
            newAnimation(obj)
        end
%%
        function obj = getNetworkData(obj,mobj)
            %could add option to select flow graph based on dispersion or
            %advections (option instead of scaling?)
            %use q for river flows and qs for drift
            caseDef = getRunParams(obj,mobj,1);
            numvar = size(obj.Data.Y,2);
            numele = length(caseDef.Element);
            rchids = unique([caseDef.Element(:).ReachID]);
            numrch = length(rchids(rchids>0));
            if numvar==numele
                graphtype = 'Network';
            elseif numvar==numrch
                graphtype = 'Reach';
            else
                obj.Data.network = [];
                return;
            end            
            [obj.Data.network,~] = Reach.getNetwork(caseDef,graphtype);
            %get node size
            obj.Data.nodesize = obj.UIset.Other;
        end  
%%
        function isn = isSnapShot(obj,option)
            %if a snap shot has been selected return true
            %option 1 - tests for a surface plot type
            %option 2 - tests fo a snap shot at point in time
            idx = find(strcmp(obj.UIset.typeList,obj.UIset.Type.String));
            idsurf = find(strcmp(obj.UIset.typeList,'surf'));
            %line plot type selected;
            isn = idx<idsurf;
            if option==2
                %snap shot selected
                issnap = isfield(obj.UIset,'Animate') && ~obj.UIset.Animate;  
                isn = isn && issnap;
            end
        end
    end
%%
%--------------------------------------------------------------------------
% Static ASM_Plots functions
%--------------------------------------------------------------------------
    methods(Static, Access=protected)
        function varorder = setVarOrder()
            %struct that holds the order of the variables for different
            %plot types        
            varnames = {'Time','Distance','Network','2D','3D','2DT','3DT',};
            %types of plot in 2,3 and 4D              
            d2 = {'Y','X'};
            d3 = {'Z','X','Y'};  
            %types of animaton in 2,3 and 4D        
            t2 = {'Y','T','X'};  
            t3 = {'Z','T','X','Y'};
            varorder = table(d2,t2,t2,d2,d3,t2,t3,'VariableNames',varnames);
        end
    end
end