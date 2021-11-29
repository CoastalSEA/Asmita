classdef ASM_PlotsUI < muiDataUI
%
%-------class help------------------------------------------------------
% NAME
%   ASM_PlotsUI.m
% PURPOSE
%   Class implements the muiDataUI to access data for use in
%   plotting
% SEE ALSO
%   muiDataUI.m, muiPLotsUI, muiPlots.m
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
% 
    properties (Transient)
        %Abstract variables for muiDataUI---------------------------        
        %names of tabs providing different data accces options
        TabOptions = {'Time','Distance','Network','2D','3D','2DT','3DT'}; 
        %selections that force a call to setVariableLists
        updateSelections = {'Case','Dataset','Variable'};
        %Additional variables for application------------------------------
        Tabs2Use         %number of tabs to include  (set in getPlotGui)     
    end  
%%  
    methods (Access=protected)
        function obj = ASM_PlotsUI(mobj)
            %initialise standard figure and menus
            guititle = 'Select Data for Plotting';
            setDataUIfigure(obj,mobj,guititle);    %initialise figure     
        end
    end
%%    
    methods (Static)
        function obj = getPlotsUI(mobj)
            %this is the function call to initialise the UI and assigning
            %to a handle of the main model UI (mobj.mUI.PlotsUI) 
            %options for selection on each tab are defined in setTabContent
            if isempty(mobj.Cases.Catalogue.CaseID)
                warndlg('No data available to plot');
                obj = [];
                return;
            elseif isa(mobj.mUI.PlotsUI,'ASM_PlotsUI')
                obj = mobj.mUI.PlotsUI;
                if isempty(obj.dataUI.Figure)
                    %initialise figure 
                    guititle = 'Select Data for Plotting';
                    setDataUIfigure(obj,mobj,guititle);    
                    setDataUItabs(obj,mobj); %add tabs 
                else
                    getdialog('Plot UI is open');
                end
            else
                obj = ASM_PlotsUI(mobj);
                idx = find(~ismember(mobj.DataUItabs.Plot,obj.TabOptions));
                if ~isempty(idx)
                    txt1 = 'Unknown plot type defined in main UI for DataUItabs.Plot';
                    for j=1:length(idx)
                        txt1 = sprintf('%s\n"%s" not defined',txt1,mobj.DataUItabs.Plot{idx(j)});
                    end
                    warndlg(txt1)
                    obj = [];
                    return
                end
                obj.Tabs2Use = mobj.DataUItabs.Plot;
                setDataUItabs(obj,mobj); %add tabs                
            end                
        end
    end
%%
%--------------------------------------------------------------------------
% Abstract methods required by muiDataUI to define tab content
%--------------------------------------------------------------------------
    methods (Access=protected) 
        function setTabContent(obj,src)
            %setup default layout options for individual tabs
            %Abstract function required by muiDataUI
            itab = find(strcmp(obj.Tabs2Use,src.Tag));
            obj.TabContent(itab) = muiDataUI.defaultTabContent;
            
            %customise the layout of each tab. Overload the default
            %template with a function for the tab specific definition
            switch src.Tag
                case 'Time'
                    setTS_tab(obj,src);
                case 'Distance'
                    setDist_tab(obj,src);
                case 'Network'
                    setNetwork_tab(obj,src);
                case '2D'
                    set2D_tab(obj,src);
                case '3D'
                    set3D_tab(obj,src);  
                case '2DT'
                    set2DT_tab(obj,src);
                case '3DT'
                    set3DT_tab(obj,src);     
            end             
        end                
%%
        function setVariableLists(obj,src,mobj)
            %initialise the variable lists or values
            %Abstract function required by muiDataUI
            %called during initialisation by muiDataUI.setDataUItabs and 
            %whenever muiDataUI.updateCaseList is called
            %Note: Scaling is handled during variable selection
            muicat = mobj.Cases.Catalogue;
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            sel_uic = S.Selections;
            caserec = sel_uic{strcmp(S.Order,'Case')}.Value;
            cobj = getCase(mobj.Cases,caserec);
            dsnames = fieldnames(cobj.Data); %**
            
            for i=1:length(sel_uic)                
                switch sel_uic{i}.Tag
                    case 'Case'                        
                        sel_uic{i}.String = muicat.CaseDescription;
                        sel_uic{i}.UserData = sel_uic{i}.Value; %used to track changes
                    case 'Dataset'
                        sel_uic{i}.String = dsnames; %**
                        sel_uic{i}.UserData = sel_uic{i}.Value; %used to track changes
                    case 'Variable' 
                        setval = sel_uic{2}.Value;  %current dataset selection
                        dataset = dsnames{setval};
                        sel_uic{i}.String = cobj.Data.(dataset).VariableDescriptions;
                        sel_uic{i}.UserData = sel_uic{i}.Value; %used to track changes
                    case 'Element'
                        setval = sel_uic{2}.Value;  %current dataset selection
                        dataset = dsnames{setval};
                        sel_uic{i}.String = cobj.Data.(dataset).Dimensions.EleName;
                        sel_uic{i}.Value = 1;   
                    case 'RunTime'
                        runtime = cobj.Data.(dataset).RowNames;
                        runtime.Format = 'dd-MMM-yyyy HH:mm:ss';
                        %shifts are need to ensure that selected value is
                        %within record so getProperty finds a record
                        startime = dateshift(runtime(1),'end','day')-minutes(1);
                        endtime = dateshift(runtime(end),'start','day')-days(1);
                        setTimeSlider(obj,src,startime,endtime);
                        sel_uic{i}.Min = deciyear(startime);
                        sel_uic{i}.Max = deciyear(endtime);
                        sel_uic{i}.Value = deciyear(startime);
                    case 'Type'
                        sel_uic{i}.String = S.Type;
                    otherwise
                        sel_uic{i}.String = 'Not set';
                end
            end        
            obj.TabContent(itab).Selections = sel_uic;
        end
%%
        function setTimeSlider(~,src,startime,endtime)
            %set the end markers and current value text for the slider
            htxt = findobj(src,'Tag','SLstart');
            htxt.String = deciyear(startime); %these are just markers
            htxt = findobj(src,'Tag','SLend');
            htxt.String = deciyear(endtime);
            htxt = findobj(src,'Tag','RunTimevalue');
            htxt.String = deciyear(startime);
        end
%%        
        function useSelection(obj,src,mobj)  
            %make use of the selection made to create a plot of selected type
            %Abstract function required by muiDataUI
            if strcmp(src.String,'Save')   %save animation to file
                saveAnimation(obj,src,mobj);
            else
                %some tabs use sub-selections of Cases and/or Variables
                %however muiDataUI already checks for subselection and adjusts
                ASM_Plots.getPlot(obj,mobj);
            end
        end 
    end
%%
%--------------------------------------------------------------------------
% Additional methods used to define tab content
%--------------------------------------------------------------------------
    methods (Access=private)
        function setTS_tab(obj,src)
            %customise the layout of the timeseries tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a Timeseries plot (line, bar, etc) select Case, Dataset, Variable and Element';
            txt2 = 'Assign to the Var button and adjust the variable range, and scaling, if required';
            txt3 = 'Select the Type of plot and then use the New, Add or Delete button to continue';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);  
            %Specification of uicontrol for each selection variable  
            %Use default lists except
            S.Titles = {'Case','Dataset','Variable','Element','Type'};  
            S.Order = {'Case','Dataset','Variable','Element','Type'};
            S.Style = {'popupmenu','popupmenu','popupmenu','popupmenu','popupmenu'};
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.1,0.14;0.3,0.14;0.5,0.14;0.7,0.14]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 1;                           %minimum number of buttons to use
            S.XYZmxvar = 0;                          %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.25,0.91,0.15];      %position for XYZ button panel
            S.XYZlabels = {'Var'};                   %default button labels
            
            %Action button specifications - use default
       
            obj.TabContent(itab) = S;                %update object
        end
%%
        function setDist_tab(obj,src)
            %customise the layout of the profiles tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.12];    %header vertical position and height
            txt1 = 'Select Case, Dataset, and Variable,adjust range and scaling, and assign to the Var button ';
            txt2 = 'For snap shot plot select Run Time. For animation change Ti button to Mv';
            txt3 = 'For X-T surface plot select surface Type, otherwise use line type (eg bar or stem)';
            txt4 = 'Use the New button to generate plot selection, or Add, Delete if snap shot plot';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s\n4 %s',txt1,txt2,txt3,txt4); 
            %Use default lists except   For an animated Distance plot, 
            S.Titles = {'Case','Dataset','Variable','Run Time','Type'};  
            S.Order = {'Case','Dataset','Variable','RunTime','Type'};
            S.Style = {'popupmenu','popupmenu','popupmenu','slider','popupmenu'};
            %distance plot can generate single line plot, animataion of 
            %line plot, or surface plot
            S.Type = {'line','bar','scatter','stem','stairs',...
                      'surf','contour','contourf','mesh','User'}; 
                  
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.06,0.08;0.22,0.08;0.38,0.08;0.54,0.08]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 1;                           %minimum number of buttons to use
            S.XYZmxvar = 0;                          %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.18,0.91,0.15];      %position for XYZ button panel
            S.XYZlabels = {'Var'};                   %default button labels
            %Action button specifications
            S.ActButNames = {'Refresh','Animate'}; %names assigned selection struct
            S.ActButText = {char(174),'Ti'};     %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.22];%positions for action buttons   
            % action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setAnimate(src,evt)'};
            % tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                           'Snap shot. Press to use animate'};  
            
            obj.TabContent(itab) = S;                %update object
        end
        
%%
        function setNetwork_tab(obj,src)
            %customise the layout of the network animation tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For an animated Network plot, select Case, Dataset, and Variable';
            txt2 = 'Assign to the Var button and adjust the variable range, and scaling, if required';
            txt3 = 'Select the Node size to use for Element display and then use the Select button to continue';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3); 
            %Use default lists except
            S.Titles = {'Case','Dataset','Variable','Node size'};   
            S.Style = {'popupmenu','popupmenu','popupmenu','slider'};
            S.Order = {'Case','Dataset','Variable','Other'};
            
            %Tab control button options
            %Use defaults
            
            %XYZ panel definition (if required)
            S.XYZnset = 1;                           %minimum number of buttons to use
            S.XYZmxvar = 0;                          %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.18,0.91,0.15];      %position for XYZ button panel
            S.XYZlabels = {'Var'};                   %default button labels
            %Action button specifications - use default
            
            obj.TabContent(itab) = S;                %update object
        end        
        
        
%%
        function set2D_tab(obj,src)
            %customise the layout of the 2D tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a cartesian plot (line, bar, etc) select variable and dimension to use for X and Y axes';
            txt2 = 'Assign variable or dimension to Var & X buttons and adjust ranges, and scaling, if required';
            txt3 = 'Use the XY button to swap the X and Y variables without needing to re-assign them';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            %Use default lists except
            
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.1,0.14;0.3,0.14;0.5,0.14;0.7,0.14]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 2;                           %minimum number of buttons to use
            S.XYZmxvar = [1,1];                      %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.25,0.91,0.25];      %position for XYZ button panel
            S.XYZlabels = {'Var','X'};               %default button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','Swap'}; %names assigned selection struct
            S.ActButText = {char(174),'XY'};     %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.41];%positions for action buttons   
            % action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setXYorder(src,evt)'};
            % tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                           'Swap from X-Y to Y-X axes'};         
            obj.TabContent(itab) = S;                %update object
        end
%%
        function set3D_tab(obj,src)
            %customise the layout of the 3D tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            
            %Header size and text
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a contour or surface plot, select a variable with at least 2 dimensions';
            txt2 = 'If the variable has more than 2 dimensions you will be prompted to select a sub-set';
            txt3 = 'Select dimensions to use for the X-Y axes (dimensions must be consistent with variable selected)';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);

            %Specification of uicontrol for each selection variable  
            S.Type = {'surf','contour','contourf','contour3','mesh','User'}; 
            
            %Tab control button options - use defaults
           
            %XYZ panel definition (if required)
            S.XYZnset = 3;                         %minimum number of buttons to use
            S.XYZmxvar = [0,1,1];                  %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.14,0.91,0.36];     %position for XYZ button panel
            S.XYZlabels = {'Var','X','Y'};         %button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh'};  %names assigned selection struct
            S.ActButText = {char(174)};   %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1];      %positions for action buttons   
            %action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)'};
            %tool tips for buttons             
            S.ActButTip = {'Refresh data list'}; 
                        
            obj.TabContent(itab) = S;              %update object
        end
%%
        function set2DT_tab(obj,src)
            %customise the layout of the 2DT tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            
            %Header size and text
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For an animation of a line, surface or volume';
            txt2 = 'Select the Variable, Time and the property to define the X-axis';
            txt3 = 'Ensure that Time and X selections are correctly matched to the dimensions of the variable';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            %Use default lists except
            
            %Tab control button options
            S.TabButText = {'Run','Save','Clear'};  %labels for tab button definition
            S.TabButPos = [0.1,0.03;0.3,0.03;0.5,0.03]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 3;                          %minimum number of buttons to use
            S.XYZmxvar = [2,1,1];                   %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.14,0.91,0.30];     %position for XYZ button panel
            S.XYZlabels = {'Var','Time','X'}; %button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','Polar'};   %names assigned selection struct
            S.ActButText = {char(174),'+'};        %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.17];    %positions for action buttons   
            %action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setPolar(src,evt)'};
            %tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                            'XY to Polar; X data in degrees'};   
     
            obj.TabContent(itab) = S;               %update object
        end   
%%
        function set3DT_tab(obj,src)
            %customise the layout of the 2DT tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            
            %Header size and text
            S.HeadPos = [0.88,0.1];    %header vertical position and height
            txt1 = 'For an animation of a line, surface or volume, select a variable with at least 2 dimensions';
            txt2 = 'Select the Time and properties to define the X and Y axes';
            txt3 = 'Ensure that Time and XY selections are correctly matched to the dimensions of the variable';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            %Specification of uicontrol for each selection variable  
            %Use default lists except for:
            S.Type = {'surf','contour','contourf','contour3','mesh','User'}; 
            
            %Tab control button options
            S.TabButText = {'Run','Save','Clear'};  %labels for tab button definition
            S.TabButPos = [0.1,0.03;0.3,0.03;0.5,0.03]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 3;                          %minimum number of buttons to use
            S.XYZmxvar = [3,1,1,1];                 %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.14,0.91,0.40];     %position for XYZ button panel
            S.XYZlabels = {'Var','Time','X','Y'}; %button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','Polar'};   %names assigned selection struct
            S.ActButText = {char(174),'+'};        %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.27];    %positions for action buttons   
            %action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setPolar(src,evt)'};
            %tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                            'XY to Polar; X data in degrees'}; 
    
            obj.TabContent(itab) = S;               %update object            
        end       
    end
end