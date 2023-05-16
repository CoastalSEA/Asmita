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
        vNumber = '3.2'
        vDate   = 'June 2023'
        modelName = 'Asmita'   
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    properties
        %GeoType defines the types of element that can be used
        GeoType = {'Channel','Tidalflat','Saltmarsh','Storage',...
                   'FloodDelta','DeltaFlat','EbbDelta',...
                   'Beachface','Shoreface','Spit'};
        %the following properties then define masks that define various
        %connectivity within the model
        LWtypes = [1,5,7,9]; %GeoTypes that vary with low water
        EXtypes = 6:10;      %GeoTypes that are external to inlet
        RCtypes = 1;         %Geotypes that define a reach - can be [1,5] or 1
        REtypes = 2:6;       %Geotypes that can belong to a reach - can be [2,3,4,6] or 2:6
        FNtypes = [2,3,4];   %GeoTypes that use fine eq concentration
                             %assume this excludes DeltaFlat????
        %NB: Rivers can only connect to 'Channel' - hard coded in River class.
        %Drift uses EXtypes to limit assignment of drift inputs.
        %Reach uses 'FloodDelta' as a key word.
        %Saltmarsh, ASMinterface and ASM_model use 'Saltmarsh' and 'Tidalflat'.        
        %Estuary uses 'Delta' to find all external delta elements ie not inc. FloodDelta. 
    end
    
    methods (Static)
        function obj = Asmita(SupressPrompts)                      
            %constructor function initialises GUI
            if nargin==1
                obj.SupressPrompts = SupressPrompts;
            end
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
            obj.ModelInputs.CSThydraulics = {'Estuary','WaterLevels'};
            obj.ModelInputs.AsmitaModel = {'Element','Estuary','WaterLevels',...
                              'RunProperties','RunConditions','EqCoeffParams'};
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT' amended in ASM_PLotUI
            obj.DataUItabs.Plot = {'Time','Distance','Network','2D','3D'}; 
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Taylor'};  
            
            modelLogo = 'Asmita_logo.jpg';  %default splash figure - edit to alternative
            initialiseUI(obj,modelLogo); %initialise menus and tabs  
            asmitaInputStruct(obj);
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
                'Rivers','Drift','Interventions','Hydraulics','Run Parameters'};                                  
            menu.Setup(1).Callback = repmat({'gcbo;'},[1,8]);
            menu.Setup(1).Separator = {'off','off','on','off','off','off','off','on'};%separator preceeds item
            % submenu for Estuary data
            menu.Setup(2).List = {'System Parameters','Water Levels', ...
                         'Dispersion','Model Constants'};
            menu.Setup(2).Callback = repmat({@obj.estuaryProps},[1,4]);

            % submenu for Element data
            menu.Setup(3).List = {'Define Elements','Element Parameters',...
                                'Load Parameters','Add','Delete'};
            menu.Setup(3).Callback = repmat({@obj.elementProps},[1,5]);
                                  
            % submenu for Saltmarsh data
            menu.Setup(4).List = {'Species Parameters','Equilibrium Marsh Depth',...
                'Biomass Distribution','Marsh-flat Animation'};
            menu.Setup(4).Callback = repmat({@obj.saltmarshProps},[1,4]);
            
            % submenu for River data
            menu.Setup(5).List = {'River Inputs','River Advection',...
                                    'River Timeseries','Delete Input'};
            menu.Setup(5).Callback = repmat({@obj.riverProps},[1,4]);
                                
            % submenu for Drift data
            menu.Setup(6).List = {'Drift Inputs','Drift Advection',...
                                    'Drift Timeseries','Delete Input'};
            menu.Setup(6).Callback = repmat({@obj.driftProps},[1,4]);
                                
            % submenu for Interventions data
            menu.Setup(7).List = {'Add or Edit','Clear','Load File','Change Sign'};
            menu.Setup(7).Callback = repmat({@obj.intervenProps},[1,4]);
            menu.Setup(7).Separator = {'off','off','off','on'};
            
            % submenu for Hydraulic data
            menu.Setup(8).List = {'Additional Parameters','Create Look-up','Summary Plot'};
            menu.Setup(8).Callback = repmat({@obj.hydraulicProps},[1,3]);

            % submenu for Run Properties data
            menu.Setup(9).List = {'Time Step','Conditions',...
                'Select Eq. Coefficients','Edit Eq. Coefficients',...
                'Edit Eq. Coefficients List'};
            menu.Setup(9).Callback = repmat({@obj.runProps},[1,5]);
            menu.Setup(9).Separator = {'off','off','on','off','on'};
            
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Check Input','Run Model','Derive Output'};
            menu.Run(1).Callback = repmat({@obj.runMenuOptions},[1,3]);
            
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics'};
            menu.Analysis(1).Callback = repmat({@obj.analysisMenuOptions},[1,2]);
            
            %% Help menu --------------------------------------------------
            menu.Help.List = {'Documentation','Manual','Theory 1','Theory 2'};
            menu.Help.Callback = repmat({@obj.Help},[1,4]);
            
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
            tabs.Cases  = {'  Cases  ',@obj.refresh};  
            tabs.Settings = {'  Settings  ',''};
            subtabs.Settings(1,:) = {' System ',@obj.InputTabSummary};
            subtabs.Settings(2,:)   = {' Elements ',@obj.InputTabSummary};
            subtabs.Settings(3,:) = {' Saltmarsh ',@obj.setAsmitaTabs};
            subtabs.Settings(4,:) = {' Interventions ',@(src,evt)Interventions.IntTabPlot(obj,src,evt)};
            subtabs.Settings(5,:)   = {' Run Parameters ',@obj.InputTabSummary};
            
            tabs.Tides   = {'  Tides  ',''};
            subtabs.Tides(1,:)   = {' Network ',@(src,evt)setgraph(obj,src,evt)};
            subtabs.Tides(2,:)   = {' Water Levels ',@obj.setAsmitaTabs};
            subtabs.Tides(3,:)   = {' Hydraulics ',@obj.setAsmitaTabs};
            subtabs.Tides(4,:)   = {' TP Discharge ',@obj.setAsmitaTabs};
            subtabs.Tides(5,:) = {' TP Network ',@(src,evt)setgraph(obj,src,evt)};
            
            tabs.Flows   = {'  Flows  ',''};
            subtabs.Flows(1,:) = {' Rivers ',@(src,evt)setgraph(obj,src,evt)};
            subtabs.Flows(2,:) = {' Drift ',@(src,evt)setgraph(obj,src,evt)};            
            subtabs.Flows(3,:) = {' Input Summary ',@(src,evt)Advection.inputSummary(obj,src,evt)};
            subtabs.Flows(4,:) = {' River Input ',@(src,evt)River.TSplot(obj,src,evt)};
            subtabs.Flows(5,:) = {' Drift Input ',@(src,evt)Drift.TSplot(obj,src,evt)};
            
            tabs.Response = {' Response ',@(src,evt)Estuary.Response(obj,src,evt)};
            tabs.Plot   = {'  Q-Plot  ',@obj.getTabData};
            tabs.Stats = {'   Stats   ',@obj.setTabAction};
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
                'Estuary','System',[0.95,0.50],{200,60},'Estuary parameters:';...
                'WaterLevels','System',[0.95,0.98],{160,80},'Hydraulic parameters:';...
                'CSThydraulics','System',[0.52,0.52],{170,100},'Additional parameters:';...
                'Element','Elements',[0.95,0.95],{180,60},'Element parameters:';...
                'Saltmarsh','Saltmarsh',[0.95,0.54],{165,120},'Saltmarsh parameters:';...
                'RunProperties','Run Parameters',[0.92,0.48],{180,60},'Run time parameters:';...
                'RunConditions','Run Parameters',[0.60,0.48],{180,60},'Run conditions (true or false):';...
                'EqCoeffParams','Run Parameters',[0.92,0.96],{180,60},'Equilibrium coefficients:'};
        end    
 %%
        function setTabAction(obj,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src)
            msg = 'No results to display';
            switch src.Tag                                    
                case 'Plot' 
                     tabPlot(cobj,src,obj);
                case 'Stats'
                    lobj = getClassObj(obj,'mUI','Stats',msg);
                    if isempty(lobj), return; end
                    tabStats(lobj,src);    
            end          
        end
%%
        function setAsmitaTabs(obj,src,evt)
            %update the SSettings tabs with tables and/or plots
            switch src.Tag
                case 'Saltmarsh'
                    InputTabSummary(obj,src,evt)
                    msgtxt = 'Saltmarsh parameters have not been defined';
                    cobj = getClassObj(obj,'Inputs','Saltmarsh',msgtxt);                    
                case 'Water Levels'
                    msgtxt = 'Water level parameters have not been defined';
                    cobj = getClassObj(obj,'Inputs','WaterLevels',msgtxt);
                case 'Hydraulics'
                    msgtxt = 'Additional hydraulic parameters have not been defined';
                    cobj = getClassObj(obj,'Inputs','CSThydraulics',msgtxt);
                case 'TP Discharge'
                    msgtxt = 'River parameters have not been defined';
                    cobj = getClassObj(obj,'Inputs','River',msgtxt);
            end
            %
            if isempty(cobj), return; end
            tabPlot(cobj,src,obj);
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
                case 'System Parameters'
                    Estuary.setInput(obj);
                    tabname = 'System';
                case 'Water Levels'
                    WaterLevels.setInput(obj);
                    tabname = 'System';
                case 'Dispersion'
                    Estuary.setDispersion(obj);
                    tabname = 'Network';
                case 'Model Constants'
                    obj.Constants = setInput(obj.Constants);
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
                case 'Element Parameters'
                    Element.setInput(obj);
                case 'Load Parameters'
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
                case 'Species Parameters'
                    Saltmarsh.setInput(obj);
                    tabUpdate(obj,'System');
                case 'Equilibrium Marsh Depth'
                    Saltmarsh.EqDepthBiomassPlot(obj);
                case 'Biomass Distribution'
                    Saltmarsh.BiomassDistributionPlot(obj);
                case 'Marsh-flat Animation'
                    Saltmarsh.MarshFlatAnimation(obj);
            end
        end
%%        
        function riverProps(obj,src,~)
            %callback functions to setup river inputs
            switch src.Text
                case 'River Inputs'
                    River.setRiver(obj);
                case 'River Advection'
                    Advection.setAdvection(obj,'River');
                case 'River Timeseries'
                    River.setRiverTS(obj);    
                case 'Delete Input'
                    River.delRiver(obj);
            end
        end
%%
        function driftProps(obj,src,~)
            %callback functions to setup drift input
            switch src.Text
                case 'Drift Inputs'
                    Drift.setDrift(obj);
                case 'Drift Advection'
                    Advection.setAdvection(obj,'Drift');
                case 'Drift Timeseries'
                    Drift.setDriftTS(obj);
                case 'Delete Input'
                    Drift.delDrift(obj);
            end
        end
%%
        function intervenProps(obj,src,~)
            %callback functions to setup interventions
            switch src.Text
                case 'Add or Edit'
                    Interventions.setInterventions(obj);
                case 'Clear'
                    Interventions.clearEleInt(obj);
                case 'Load File'
                    Interventions.loadInterventions(obj);
                case 'Change Sign'
                    Interventions.setIntChSign(obj)
                    %what about add and delete?????
            end
        end
%%
        function hydraulicProps(obj,src,~)
            %callback functions to setup hydraulics
            msgtxt = ('Hydraulic properties have not been defined');
            switch src.Text
                case 'Additional Parameters'
                    CSThydraulics.setInput(obj);
                    tabname = 'System';
                    tabUpdate(obj,tabname);
                case 'Create Look-up'                    
                    cobj = getClassObj(obj,'Inputs','CSThydraulics',msgtxt);
                    if isempty(cobj), return; end
                    runModel(cobj,obj);
                case 'Summary Plot'
                    cobj = getClassObj(obj,'Inputs','CSThydraulics',msgtxt);
                    if isempty(cobj), return; end
                    tabPlot(cobj);
            end
        end
%%
        function runProps(obj,src,~)
            %set-up various types of run property
            tabname = 'Run Parameters';
            switch src.Text
                case 'Time Step'
                    RunProperties.setInput(obj);                    
                case 'Conditions'
                    RunConditions.setInput(obj);   
                case 'Select Eq. Coefficients'
                    EqCoeffParams.setInput(obj);
                case 'Edit Eq. Coefficients'  
                    EqCoeffParams.editParams(obj);
                case 'Edit Eq. Coefficients List'
                    EqCoeffParams.editList(obj);
                    tabname = [];
            end
            %                
            tabUpdate(obj,tabname);
        end
        %% Run menu -------------------------------------------------------
        function runMenuOptions(obj,src,~)
            %callback functions to run model
            switch src.Text        
                case 'Check Input'
                    [message,ok] = AsmitaModel.CheckInput(obj);
                    msgbox(message,'Check Inputs');
                    if ok<1
                        warndlg('Not enough components defined to run model');
                    end
                case 'Run Model'                        
                    AsmitaModel.runModel(obj); 
                case 'Derive Output'
                    obj.mUI.ManipUI = muiManipUI.getManipUI(obj);
            end            
        end               
            
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = ASM_PlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,src,~)
            switch src.Text
                case 'Documentation'
                    doc asmita   
                case 'Manual'
                    asm_open_manual;
                case 'Theory 1'
                    asm_coe_part1;
                case 'Theory 2'
                    asm_coe_part2;
            end
        end    
    end
%%
    methods
        function loadAsmitaModel(obj,path,file)
            %load a model (used in unit testing)
            obj.Info.PathName = path;
            obj.Info.FileName = file;
            loadModel(obj);
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
%%
        function asmitaInputStruct(obj)
            %define struct to be used as handles for input data classes
            obj.Inputs = struct('WaterLevels',[],'Estuary',[],'Element',[],...
             'RunProperties',[],'RunConditions',[],'EqCoeffParams',[],...
             'River',[],'Drift',[],'Saltmarsh',[],'Interventions',[],...
             'Reach',[],'Advection',[],'CSThydraulics',[]);
        end
    end
end    
    
    
    
    
    
    
    
    
    
    