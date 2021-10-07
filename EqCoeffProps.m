classdef EqCoeffProps < muiPropertyUI                
%
%-------class help---------------------------------------------------------
% NAME
%   EqCoeffProps.m
% PURPOSE
%   Class for setting the selected equilibrium coefficients
% USAGE
%   obj = EqCoeffProps.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      

%cannot use muiProperties because array is n-elements x 3 properties
%Strip eqcoeffs out of old RunProps


    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Alpha','Beta','EqType'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        %list of options used in UserEqPrismCoeffs
        UserEqCoeffOptions ={'Default','Generic','Venice','Amelander','Pagham'};
        UserEqCoeffSelection %selection used in UserEqPrismCoeffs (list defined in AsmitaModel)
        SupressPrompts = false %flag for unit testing to supress user promts
    end
    
    properties
        alpha                %scale coefficients for equilibrium prism relationship
        beta                 %shape coefficients for equilibrium prism relationship
        eqtype               %use tidal prism or tidal range to determine equilibrium 
    end    
    

%%   
    methods (Access={?muiPropertyUI,?muiModelUI})
        function obj = EqCoeffProps(mobj)          
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'EqCoeffProps'; 
            obj = EqCoeffProps.setClassObj(mobj,classname);
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                selections = obj.UserEqCoeffOptions;
                h_dlg = listdlg('PromptString','Select a coefficient set:',...
                    'SelectionMode','single','ListString',selections);
                if isempty(h_dlg), return; end
                obj.UserEqCoeffSelection = selections{h_dlg};
                [obj.alpha,obj.beta,obj.eqtype] = UserPrismCoeffs(selections{h_dlg});
            end
            mobj.Inputs.(classname) = obj;
        end 
%%
        function obj = editList(mobj)            
            %create table to input or edit equlibrium coefficients list
            classname = 'EqCoeffProps'; 
            obj = EqCoeffProps.setClassObj(mobj,classname);
            userdata = obj.UserEqCoeffOptions';
            oldtable = table(userdata,'VariableNames',{'EqCoeff_Options'});
                                    
            title = sprintf('Edit Equilibrium List'); 
            header = 'Select row to edit or use Add button for new definition';  
            but.Text = {'Save','Add','Cancel'}; %labels for tab button definition
            newtable = tablefigureUI(title,header,oldtable,true,but);
            if isempty(newtable), return; end  %user cancelled  
            obj.UserEqCoeffOptions =  newtable.EqCoeff_Options;
            mobj.Inputs.(classname) = obj;
        end           
    end
%%        
        %add other functions to operate on properties as required  
        
    methods
        function displayProperties(obj,src)
            %table for equilibirum coefficient selection alpha and beta.
            %called by muiModelUI.InputTabSummary and overloads super class
            %muiPropertyUI method.
            if isempty(obj.UserEqCoeffSelection)
                obj.UserEqCoeffSelection = 'Not set';
            end
            if isstruct(obj.alpha)                
%                 rownames = unique(eletype);
                rownames = fieldnames(obj.alpha);
                colnames = {'Type','Alpha','Beta','EqType'};
                ntypes = length(rownames);
                userdata = cell(ntypes,3);
                for i=1:ntypes
                    userdata{i,1} = rownames{i};
                    userdata{i,2} = obj.alpha.(rownames{i});
                    userdata{i,3} = obj.beta.(rownames{i});
                    userdata{i,4} = obj.eqtype.(rownames{i});
                end
                
                Table = uitable('Parent',src, ...
                    'ColumnName', colnames, ...
                    'RowName', [], ....
                    'ColumnWidth', {80 60 55 45}, ...
                    'Data',userdata, ...
                    'Units','normalized');
                Table.Position(3:4)=Table.Extent(3:4);
                Table.Position(2)=0.9-Table.Extent(4);
                Table.Position(1)=0.95-Table.Extent(3);
            end
            helptxt = sprintf('Equilibrium Coefficients using: %s',obj.UserEqCoeffSelection);
            uicontrol('Parent',src,...
                    'Style','text','String', helptxt,...                    
                    'HorizontalAlignment', 'left',...
                    'Units','normalized', 'Position', [0.55 0.9 0.4 0.05],...
                    'Tag','RPtext');
        end  
     
    end
%%
    methods (Static, Access=private)
        function obj = setClassObj(mobj,classname)
            %check if class exists and if not call constructor                
            if isfield(mobj.Inputs,classname) && ...
                            isa(mobj.Inputs.(classname),classname)
                obj = mobj.Inputs.(classname);  
            else
                obj = EqCoeffProps(mobj);    
            end
        end
    end
end