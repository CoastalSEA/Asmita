classdef Asmita < muiModelUI                      
%
%-------class help---------------------------------------------------------
% NAME
%   Asmita.m
% PURPOSE
%   Main GUI for a generic model interface, which implements the 
%   muiModelUI abstract class to define main menus.
% SEE ALSO
%   Abstract class muiModelUI.m and tools provided in muitoolbox
%
% Author: Ian Townend
% CoastalSEA (c) Apr 2021
%--------------------------------------------------------------------------
% 
    properties  (Access = protected)
        %implement properties defined as Abstract in muiModelUI
        vNumber = '3.0'
        vDate   = 'May 2021'
        modelName = 'Asmita'   
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    properties
        %GeoType defines the types of element that have been implemented
        GeoType = {'Delta','Channel','Tidalflat','Saltmarsh','Storage',...
                   'EbbDelta','FloodDelta','DeltaFlat',...
                   'Beachface','Shoreface','Spit'};
    end
    
    methods (Static)
        function obj = Asmita                      
            %constructor function initialises GUI
            obj = setMUI(obj);             
        end
    end
%% ------------------------------------------------------------------------
% Definition of GUI Settings
%--------------------------------------------------------------------------  
    methods (Access = protected)
        function obj = setMUI(obj)
            %initialise standard figure and menus    
            %classes required to run model             
            %format:                                         % << Edit to model and input parameters classnames 
            %obj.ModelInputs.<model classname> = {'Param_class1',Param_class2',etc}
            obj.ModelInputs.Model_template = {'ParamInput_template'};
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'2D','3D','4D','2DT','3DT','4DT'};  
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Timeseries','Taylor','Intervals'};  
            
            modelLogo = 'Asmita_logo.jpg';  %default splash figure - edit to alternative
            initialiseUI(obj,modelLogo); %initialise menus and tabs                  
        end    
        
%% ------------------------------------------------------------------------
% Definition of Menu Settings
%--------------------------------------------------------------------------
        function menu = setMenus(obj)
            %define top level menu items and any submenus
            %MenuLabels can any text but should avoid these case-sensitive 
            %reserved words: "default", "remove", and "factory". If label 
            %is not a valid Matlab field name this the struct entry
            %is modified to a valid name (eg removes space if two words).
            %The 'gcbo:' Callback text triggers an additional level in the 
            %menu. Main menu labels are defined in sequential order and 
            %submenus in order following each brach to the lowest level 
            %before defining the next branch.         
                                                              % << Edit menu to suit model 
            MenuLabels = {'File','Tools','Project','Setup','Run',...
                                                        'Analysis','Help'};
            menu = menuStruct(obj,MenuLabels);  %create empty menu struct
            %
            %% File menu --------------------------------------------------
             %list as per muiModelUI.fileMenuOptions
            menu.File.List = {'New','Open','Save','Save as','Exit'};
            menu.File.Callback = repmat({@obj.fileMenuOptions},[1,5]);
            
            %% Tools menu -------------------------------------------------
            %list as per muiModelUI.toolsMenuOptions
            menu.Tools(1).List = {'Refresh','Clear all'};
            menu.Tools(1).Callback = {@obj.refresh, 'gcbo;'};  
            
            % submenu for 'Clear all'
            menu.Tools(2).List = {'Model','Figures','Cases'};
            menu.Tools(2).Callback = repmat({@obj.toolsMenuOptions},[1,3]);

            %% Project menu -----------------------------------------------
            menu.Project(1).List = {'Project Info','Cases','Export/Import'};
            menu.Project(1).Callback = {@obj.editProjectInfo,'gcbo;','gcbo;'};
            
            %list as per muiModelUI.projectMenuOptions
            % submenu for Scenarios
            menu.Project(2).List = {'Edit Description','Edit Data Set',...
                                    'Save Data Set','Delete Case','Reload Case',...
                                    'View Case Settings'};                                               
            menu.Project(2).Callback = repmat({@obj.projectMenuOptions},[1,6]);
            
            % submenu for 'Export/Import'                                          
            menu.Project(3).List = {'Export Case','Import Case'};
            menu.Project(3).Callback = repmat({@obj.projectMenuOptions},[1,2]);
            
            %% Setup menu -------------------------------------------------                     
            menu.Setup(1).List = {'Estuary','Elements','Saltmarsh',...
                'Rivers','Drift','Interventions','Run Parameters'};                                  
            menu.Setup(1).Callback = repmat({'gcbo;'},[1,7]);
            menu.Setup(1).Separator = {'off','off','on','off','off','off','on'};%separator preceeds item
            % submenu for Estuary data
            menu.Setup(2).List = {'System Properties','Water Levels', ...
                         'Dispersion','Hydraulic model','Model Constants'};
            menu.Setup(2).Callback = repmat({@obj.estuaryProps},[1,5]);

            % submenu for Element data
            menu.Setup(3).List = {'Define Elements','Element Properties',...
                                'Load Properties','Add','Delete'};
            menu.Setup(3).Callback = repmat({@obj.elementProps},[1,5]);
                                  
            % submenu for Saltmarsh data
            menu.Setup(4).List = {'Species parameters','Equilibrium marsh depth',...
                'Biomass distribution','Marsh-flat animation'};
            menu.Setup(4).Callback = repmat({@obj.saltmarshProps},[1,4]);
            
            % submenu for River data
            menu.Setup(5).List = {'River inputs','River advection',...
                                    'River timeseries','Delete input'};
            menu.Setup(5).Callback = repmat({@obj.riverProps},[1,4]);
                                
            % submenu for Drift data
            menu.Setup(6).List = {'Drift inputs','Drift advection',...
                                    'Drift timeseries','Delete input'};
            menu.Setup(6).Callback = repmat({@obj.driftProps},[1,4]);
                                
            % submenu for Interventions data
            menu.Setup(7).List = {'Add','Load'};
            menu.Setup(7).Callback = {@obj.interProps,@obj.interProps};
            
            % submenu for Run Properties data
            menu.Setup(8).List = {'Time step','Conditions',...
                'Equilibrium Coefficients','Edit Eq. Coefficients list'};
            menu.Setup(8).Callback = repmat({@obj.runProps},[1,4]);
            
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Run Model','Derive Output'};
            menu.Run(1).Callback = repmat({@obj.runMenuOptions},[1,2]);
            
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics'};
            menu.Analysis(1).Callback = repmat({@obj.analysisMenuOptions},[1,2]);
            
            %% Help menu --------------------------------------------------
            menu.Help(1).Callback = {@obj.Help}; %make model specific?
            
        end
        
%% ------------------------------------------------------------------------
% Definition of Tab Settings
%--------------------------------------------------------------------------
        function [tabs,subtabs] = setTabs(obj)
            %define main tabs and any subtabs required. struct field is 
            %used to set the uitab Tag (prefixed with sub for subtabs). 
            %Order of assignment to struct determines order of tabs in figure.
            %format for tabs: 
            %    tabs.<tagname> = {<tab label>,<callback function>};
            %format for subtabs: 
            %    subtabs.<tagname>(i,:) = {<subtab label>,<callback function>};
            %where <tagname> is the struct fieldname for the top level tab.
            tabs.Cases  = {'  Cases  ',@obj.refresh};        % << Edit tabs to suit model 
            tabs.System = {' System ',@obj.InputTabSummary};
            tabs.Elements   = {' Elements ',@obj.InputTabSummary};
            tabs.Interven = {' Interven ',@(src,evdat)Interventions.IntSummary(obj,src,evdat)};
            tabs.RunProps   = {' Run Props ',@obj.InputTabSummary};
            tabs.Network   = {' Network ',@(src,evt)Estuary.Network(obj,src,evt)};
            tabs.Flows   = {'  Flows  ','gcbo;'};
            subtabs.Flows(1,:) = {' Rivers ',@(src,evt)Estuary.Network(obj,src,evt)};
            subtabs.Flows(2,:) = {' Drift ',@(src,evt)Estuary.Network(obj,src,evt)};
            subtabs.Flows(3,:) = {'Tidal Pumping',@(src,evt)Estuary.Network(obj,src,evt)};
            subtabs.Flows(4,:) = {'Input Summary',@(src,evt)Advection.inputSummary(obj,src,evt)};
            subtabs.Flows(5,:) = {'River Input',@(src,evt)River.TSplot(obj,src,evt)};
            subtabs.Flows(6,:) = {'Drift Input',@(src,evt)Drift.TSplot(obj,src,evt)};
            tabs.Response = {' Response ',@(src,evdat)Estuary.Response(obj,src,evdat)};            
        end
       
%%
        function props = setTabProperties(~)
            %define the tab and position to display class data tables
            %props format: {class name, tab tag name, position, ...
            %               column width, table title}
            % position and column widths vary with number of parameters
            % (rows) and width of input text and values. Inidcative
            % positions:  top left [0.95,0.48];    top right [0.95,0.97]
            %         bottom left [0.45, 0.48]; bottom rigth [0.45,0.97]                                                    
            props = {...                                     
                'Estuary','System',[0.95,0.50],{200,60},'Estuary properties:';...
                'WaterLevels','System',[0.95,0.98],{180,60},'Hydraulic properties:';...
                'Saltmarsh','System',[0.55,0.55],{170,100},'Saltmarsh properties:';...
                'Element','Elements',[0.95,0.95],{180,60},'Element properties:';...
                'RunProps','RunProps',[0.95,0.48],{180,60},'Run time properties:';...
                'RunConditions','RunProps',[0.55,0.48],{180,60},'Run conditions (true or false):';...
                'EqCoeffProps','RunProps',[0.95,0.95],{180,60},'Equilibrium coefficients:'};
        end    
 %%
        function setTabAction(~,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src)
            switch src.Tag                                    % << Edit match tab requirements
                case 'Plot' 
                     tabPlot(cobj,src);
                case 'Stats'
                     tabStats(cobj,src);    
            end
        end      
%% ------------------------------------------------------------------------
% Callback functions used by menus and tabs
%-------------------------------------------------------------------------- 
        %% File menu ------------------------------------------------------
        %use default menu functions defined in muiModelUI
            
        %% Tools menu -----------------------------------------------------
        %use default menu functions defined in muiModelUI
                
        %% Project menu ---------------------------------------------------
        %use default menu functions defined in muiModelUI           

        %% Setup menu -----------------------------------------------------     
        function estuaryProps(obj,src,~)
            %call functions to setup estuary definition
            tabname = [];
            switch src.Text
                case 'System Properties'
                    Estuary.setInput(obj);
                    tabname = 'System';
                case 'Water Levels'
                    WaterLevels.setInput(obj);
                    tabname = 'System';
                case 'Dispersion'
                    Estuary.setDispersion(obj);
                    tabname = 'Network';
                case 'Hydraulic model'
                case 'Model Constants'
                    obj.Constants = editProperties(obj.Constants);
            end
            %
            tabUpdate(obj,tabname);
        end
%%
        function elementProps(obj,src,~)
            %callback functions to setup elements
            switch src.Text
                case 'Define Elements'
                    Element.defineElements(obj);
                case 'Element Properties'
                    Element.setInput(obj);
                case 'Load Properties'
                    Element.loadElement(obj);
                case 'Add'
                    Element.addElement(obj);
                case 'Delete'
                    Element.delElement(obj);
            end
        end   
%% 
        function saltmarshProps(obj,src,~)
            %callback functions to setup saltmarsh
            switch src.Text
                case 'Species parameters'
                    Saltmarsh.setInput(obj);
                    tabUpdate(obj,'System');
                case 'Equilibrium marsh depth'
                    Saltmarsh.EqDepthBiomass(obj);
                case 'Biomass distribution'
                    Saltmarsh.BiomassDistribution(obj);
                case 'Marsh-flat animation'
                    Saltmarsh.MarshFlatAnimation(obj);
            end
        end
%%        
        function riverProps(obj,src,~)
            %callback functions to setup river inputs
            switch src.Text
                case 'River inputs'
                    River.setRiver(obj);
                case 'River advection'
                    Advection.setAdvection(obj,'River');
                case 'River timeseries'
                    River.setRiverTS(obj);    
                case 'Delete input'
                    River.delRiver(obj);
            end
        end
%%
        function driftProps(obj,src,~)
            %callback functions to setup drift input
            switch src.Text
                case 'Drift inputs'
                    Drift.setDrift(obj);
                case 'Drift advection'
                    Advection.setAdvection(obj,'Drift');
                case 'Drift timeseries'
                    Drift.setDriftTS(obj);
                case 'Delete input'
                    Drift.delDrift(obj);
            end
        end
%%
        function intervenProps(obj,src,~)
            %callback functions to setup interventions
            switch src.Text
                case 'Add'
                    Interventions.setInterventions(obj);
                case 'Load'
                    Interventions.loadInterventions(obj);
            end
        end
%%
        function runProps(obj,src,~)
            %set-up various types of run property
            tabname = 'RunProps';
            switch src.Text
                case 'Time step'
                    RunProps.setInput(obj);                    
                case 'Conditions'
                    RunConditions.setInput(obj);   
                case 'Equilibrium Coefficients'
                    cobj = setInputClassObj(mobj,'EqCoeffProps');
                    setInput(cobj);
                case 'Edit Eq. Coefficients list'                    
                    EqCoeffProps.editList(obj);
                    tabname = [];
            end
            %
            tabUpdate(obj,tabname);
        end
        %% Run menu -------------------------------------------------------
        function runMenuOptions(obj,src,~)
            %callback functions to run model
            switch src.Text                   
                case 'Run Model'                             % << Edit to call Model class
                    Model_template.runModel(obj); 
                case 'Derive Output'
                    obj.mUI.Manip = muiManipUI.getManipUI(obj);
            end            
        end               
            
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.Plots = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.Stats = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,~,~)
            doc UseUI_template                               % << Edit to documentation name if available
        end    
    end
%%
    methods (Access=private)
        function tabUpdate(obj,tabname)
            %update tab used for properties if required
            if ~isempty(tabname)
                tabsrc = findobj(obj.mUI.Tabs,'Tag',tabname);
                InputTabSummary(obj,tabsrc);
            end 
        end
    end
end    
    
    
    
    
    
    
    
    
    
    