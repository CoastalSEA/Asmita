classdef EqCoeffParams < muiPropertyUI                
%
%-------class help---------------------------------------------------------
% NAME
%   EqCoeffParams.m
% PURPOSE
%   Class for setting the selected equilibrium coefficients
% USAGE
%   obj = EqCoeffParams.setInput(mobj); %mobj is a handle to Main UI
% NOTES
%   table displays element types x 3 properties
% SEE ALSO
%   inherits muiPropertyUI uses definitions in userprismcoeffs function
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Alpha','Beta','EqType'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        %list of options used in UserEqPrismCoeffs
        UserEqCoeffOptions ={'Default','Generic-linear','Generic-power',...
                                            'Venice','Amelander','Pagham'};
        UserEqCoeffSelection %selection used in UserEqPrismCoeffs (list defined in AsmitaModel)
        SupressPrompts = false %flag for unit testing to supress user prompts
    end
    
    properties
        alpha     %scale coefficients for equilibrium prism relationship
        beta      %shape coefficients for equilibrium prism relationship
        eqtype    %use tidal prism or tidal range to determine equilibrium 
    end    
    

%%   
    methods (Access={?muiPropertyUI,?muiModelUI})
        function obj = EqCoeffParams(mobj)          
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
            classname = 'EqCoeffParams'; 
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = EqCoeffParams(mobj);  
            end
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                selections = obj.UserEqCoeffOptions;
                h_dlg = listdlg('PromptString','Select a coefficient set:',...
                    'SelectionMode','single','ListString',selections);
                if isempty(h_dlg), return; end
                obj.UserEqCoeffSelection = selections{h_dlg};
                [obj.alpha,obj.beta,obj.eqtype] = userprismcoeffs(selections{h_dlg});
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end 
%%
        function editParams(mobj)
            %edit the selected parameters by element type
            classname = 'EqCoeffParams'; 
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = EqCoeffParams(mobj);  
            end
            %
            if isempty(obj.alpha)
                obj.UserEqCoeffSelection = obj.UserEqCoeffOptions{1};
                [obj.alpha,obj.beta,obj.eqtype] = userprismcoeffs(obj.UserEqCoeffSelection);
            end
            %find code to select element by type
            typeList = fieldnames(obj.alpha);
            ok = 1;
            while ok>0
                [idx, ok] = listdlg('Name','EqCoeffs',...
                'PromptString','Select element type','SelectionMode','single',...
                'ListString',typeList);
                if ok==0, continue; end   %Use cancel to quit loop
                
                %edit values of selected Element Type
                prompt = {'Alpha','Beta','EqType'};
                title = 'EqCoeffs';
                numlines = 1;
                defaultvalues = {num2str(obj.alpha.(typeList{idx})),...
                                 num2str(obj.beta.(typeList{idx})),...
                                 num2str(obj.eqtype.(typeList{idx}))};                                                 
                %use updated properties to call inpudlg and return new values
                useInp=inputdlg(prompt,title,numlines,defaultvalues);
                if isempty(useInp), continue; end %user cancelled
                obj.alpha.(typeList{idx}) = str2double(useInp{1});
                obj.beta.(typeList{idx}) = str2double(useInp{2});
                obj.eqtype.(typeList{idx}) = str2double(useInp{3});
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end
%%
        function obj = editList(mobj)            
            %create table to input or edit equlibrium coefficients list
            classname = 'EqCoeffParams'; 
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = EqCoeffParams(mobj);  
            end
            userdata = obj.UserEqCoeffOptions';
            oldtable = table(userdata,'VariableNames',{'EqCoeff_Options'});
                       
            title = sprintf('Edit Equilibrium List'); 
            header = 'Select row to edit or use Add button for new definition';  
            but.Text = {'Save','Add','Cancel'}; %labels for tab button definition
            newtable = tablefigureUI(title,header,oldtable,true,but);
            if isempty(newtable), return; end  %user cancelled  
            obj.UserEqCoeffOptions =  newtable.EqCoeff_Options';
            setClassObj(mobj,'Inputs',classname,obj);
        end           
    end
%%        
        %add other functions to operate on properties as required  
        
    methods
        function displayProperties(obj,src)
            %table for equilibrium coefficient selection alpha and beta.
            %called by muiModelUI.InputTabSummary and overloads super class
            %muiPropertyUI method.
            if isempty(obj.UserEqCoeffSelection)
                obj.UserEqCoeffSelection = 'Not set';
            end
            
            tabpos = obj.TabDisplay.Position;
            if isstruct(obj.alpha)                
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

                fmtdata = @(x) sprintf('% .3e',x);
                userdata(:,2) = cellfun(fmtdata, userdata(:,2), 'UniformOutput',0);
                
                Table = uitable('Parent',src, ...
                    'ColumnName', colnames, ...
                    'RowName', [], ....
                    'ColumnWidth', {80 60 55 45}, ...
                    'Data',userdata, ...
                    'Units','normalized');
                Table.Position(3:4) = Table.Extent(3:4);
                Table.Position(2) = tabpos(1)-Table.Extent(4);
                Table.Position(1) = tabpos(2)-Table.Extent(3);
            end

            htxt = findobj(src,'Tag','RPtext');
            if isempty(htxt)
                helptxt = sprintf('Equilibrium Coefficients using: %s',obj.UserEqCoeffSelection);
                uicontrol('Parent',src,...
                    'Style','text','String', helptxt,...                    
                    'HorizontalAlignment', 'left',...
                    'Units','normalized', 'Position',[Table.Position(1),tabpos(1),0.4,0.04],...
                    'Tag','RPtext');
            end
        end       
    end
end