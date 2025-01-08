classdef CSThydraulics < muiPropertyUI                
%
%-------class help---------------------------------------------------------
% NAME
%   CSThydraulics.m
% PURPOSE
%   Class for input parameters and running the CSTmodel to be use as input 
%   to the definition of Reaches in Asmita (Setup>Hydraulics menu)
% USAGE
%   obj = CSThydraulics.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Estuary length (m)',...
                          'Width at mouth (m)',...
                          'Area at mouth (m^2)',...
                          'River width (m)',...
                          'River Area (m^2)',...
                          'Distance from mouth to estuary/river switch (m)',...
                          'Manning friction coefficient [mouth switch head]',...
                          'Storage width ratio [mouth switch head]',...
                          'River discharge range [Q1 Q2 ...  Qn]'}
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        EstuaryLength  %estuary length (m)
        MouthWidth     %width at mouth (m)
        MouthCSA       %area at mouth (m^2)    
        RiverWidth     %upstream river width (m) 
        RiverCSA       %upstream river cross-sectional area (m^2)
        xTideRiver     %distance from mouth to estuary/river switch
        Manning        %Manning friction coefficient [mouth switch head]
        StorageRatio   %storage width ratio [mouth switch head]
        Qrange         %vector of input river discharges (m^3/s)
        
    end
    
    properties (Hidden)
        CSTmodel         %dstable from CSTrunmodel with following fields
        % MeanTideLevel  %mean water suface elevation along estuary
        % TidalElevAmp   %elevation amplitude along estuary
        % TidalVelAmp    %tidal velocity amplitude along estuary
        % RiverVel       %river flow velocity along estuary
        % HydDepth       %along estuary hydraulic depth 
    end   
    
    properties (Transient)
        msgtxt = 'To include hydraulics in Asmita the CSTmodel App must be installed';
    end

%%   
    methods (Access=protected)
        function obj = CSThydraulics(mobj)
            %constructor code:
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
            %check that CSTmodel is available as an App
            ok = initialise_mui_app('CSTmodel',obj.msgtxt,'CSTfunctions');
            if ok<1, return; end
        end
    end
    %%
    methods (Static)
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'CSThydraulics';
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CSThydraulics(mobj);
            end

            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end
    end
    %%
    methods
        function runModel(obj,mobj)
            %compile input data and run model
            ok = initialise_mui_app('CSTmodel',obj.msgtxt,'CSTfunctions');
            if ok<1, return; end
            
            if isempty(obj.MouthWidth)  %check first property has been set
                warndlg('Hydraulic properties have not been defined')
                return;
            end
            dsp = modelDSproperties(obj);

            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)
                warndlg('Use Setup to define Estuary and/or Water Level parameters');
                return;
            end
            %--------------------------------------------------------------------------
            % Model code
            %--------------------------------------------------------------------------
            %input parameters for model
            [inp,rnp] = getModelParameters(obj,mobj);
            est = [];  %observed values of estuary form so can be empty
            if isempty(obj.Qrange) || all(obj.Qrange==0)
                obj.Qrange = obj.RiverDischarge;
            end

            if length(inp.AreaELength)>1   %linear variation of convergence
                LA = @(Q) inp.AreaELength(1)+inp.AreaELength(2)*Q;
            else
                LA = @(Q) inp.AreaELength; %constant convergence length
            end

            nrow = length(obj.Qrange);
            resX{nrow,5} = [];
            for i=1:nrow
                inp.RiverDischarge = obj.Qrange(i);
                inp.AreaELength = LA(inp.RiverDischarge);
                try
                    [res,xdim] = cst_model(inp,rnp,est);
                    resX(i,:) = [res.X(1,[1,2,4,5]),res.F(1,2)];
                catch ME
                    %remove the waitbar if program did not complete
                    model_catch(ME,'cst_model','Qr',inp.RiverDischarge); 
                end
            end
            resXQ = cell(1,5);
            for col = 1:5
               resXQ{col} = vertcat(resX{:,col});
            end
            %now assign results to object properties
            %--------------------------------------------------------------------------
            % Assign model output to a dstable using the defined dsproperties meta-data
            %--------------------------------------------------------------------------
            dst = dstable(resXQ{:},'RowNames',inp.Qrange','DSproperties',dsp);
            dst.Dimensions.X = xdim{1};     %grid x-coordinate
            %--------------------------------------------------------------------------
            % Save results
            %--------------------------------------------------------------------------
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            %save results
            obj.CSTmodel = dst;
            setClassObj(mobj,'Inputs','CSThydraulics',obj);
            getdialog('Run complete');
        end
%%
        function tabPlot(obj,src,~)
            %create a summary of the hydraulic model results
            if isempty(obj.CSTmodel)
                warndlg('No model results to display')
                return;
            end
            
            if nargin<2
                src = figure('Name','Hydraulic plot','Tag','PlotFig');
            else
                ax = findobj(src,'Tag','PlotFigAxes');
                delete(ax)
            end
            ax = axes('Parent',src,'Tag','PlotFigAxes');
            ax.Position = [0.16,0.18,0.65,0.75]; %make space for slider bar
            setYaxisLimits(obj,ax);  %wakes dynamic properties if necessary
            
            Q = obj.Qrange;
            setSlideControl(obj,src,Q(1),Q(end));
            
            cstPlot(obj,ax,Q(1))
        end
%%
        function wl = assignCSTproperties(obj,rchobj,famp,omsl,Q)
            %assign the CSTmodel properties to Asmita model elements
            % famp - factor for amplitude at mouth at time t relative to t=0
            % omsl - offset for changes in msl at mouth (eg for slr) at time t relative to t=0
            % Q - river discharge at time t
            reachLen = [rchobj(:).CumulativeLength];
            %assume constant amplitude at mouth to save running the model
            %at each time step (could replace this with a call to cst_model)
            cst = interpCSTproperties(obj,Q);
            mwl = interp1(cst.x,cst.z,reachLen,'linear');
            amp = interp1(cst.x,cst.a,reachLen,'linear');
            rchmwl = mwl(2:end)-diff(mwl)/2;
            rchamp = amp(2:end)-diff(amp)/2;
            %assume changes in amp and msl apply throughout
            wl.mwl = [mwl(1),rchmwl]+omsl;
            wl.amp = [amp(1),rchamp]*famp;
            %resultant high and low water levels
            wl.hwl = wl.mwl+wl.amp;
            wl.lwl = wl.mwl-wl.amp;
        end
    end
%%    
    methods (Access = private)
        function cst = interpCSTproperties(obj,Q)
            %interpolate the CSTmodel output for given value of river discharge
            dst = obj.CSTmodel;
            Qrows = dst.RowNames;
            cst.x = dst.Dimensions.X; 
            if ~isprop(dst,'MeanTideLevel') 
                %check if dynamic properties are active
                dst = activatedynamicprops(dst);
            end
            cst.z = interp1(Qrows,dst.MeanTideLevel,Q,'linear'); %mean tide level
            cst.a = interp1(Qrows,dst.TidalElevAmp,Q,'linear');  %tidal amplitude
			cst.U = interp1(Qrows,dst.TidalVelAmp,Q,'linear');   %tidal velocity amplitude
			cst.v = interp1(Qrows,dst.RiverVel,Q,'linear');      %river velocity 
			cst.d = interp1(Qrows,dst.HydDepth,Q,'linear');      %hydraulic depth
        end
%%
        function cstPlot(obj,ax,Q)
            %plot the along channel variables from the CSTmodel (as per tab)
            cst = interpCSTproperties(obj,Q);
            green = mcolor('green');
            orange = mcolor('orange');
			yyaxis(ax,'left')
            cla                                  %clear any existing plot lines
            plot(ax,cst.x,cst.z,'-r','LineWidth',1.0);           %plot time v elevation
            hold on
            plot(ax,cst.x,(cst.z+cst.a),'-.b','LineWidth',0.8)   %plot high water level
            plot(ax,cst.x,(cst.z-cst.a),'-.b','LineWidth',0.8)   %plot low water level
			plot(ax,cst.x,(cst.z-cst.d),'-k','LineWidth',0.6);   %hydraulic depth below mean tide level
            ylabel('Elevation (mOD)'); 
			yyaxis(ax,'right')
            cla                                  %clear any existing plot lines
			plot(ax,cst.x,cst.U,'--','Color',orange,'LineWidth',0.6)%plot tidal velocity
			plot(ax,cst.x,cst.v,'--','Color',green,'LineWidth',0.6) %plot river velocity
            hold off
            xlabel('Distance from mouth (m)'); 
            ylabel('Velocity (m/s)'); 
			legend('MTL','HWL','LWL','Hydraulic depth',...
                'Tidal velocity','River velocity','Location','best');			
            title('Along channel variation');
        end

%%
        function hm = setSlideControl(obj,hfig,qmin,qmax)
            %intialise slider to set different Q values   
            invar = struct('sval',[],'smin',[],'smax',[],'size', [],...
                           'callback','','userdata',[],'position',[],...
                           'stxext','','butxt','','butcback','');            
            invar.sval = qmin;     %initial value for slider 
            invar.smin = qmin;     %minimum slider value
            invar.smax = qmax;     %maximum slider value
            invar.callback = @(src,evt)updateCSTplot(obj,src,evt); %callback function for slider to use
            invar.position = [0.15,0.005,0.45,0.04]; %position of slider
            invar.stext = 'River discharge = ';   %text to display with slider value, if included          
            hm = setfigslider(hfig,invar);   
        end   
%%
        function updateCSTplot(obj,src,~)
            %use the updated slider value to adjust the CST plot
            stxt = findobj(src.Parent,'Tag','figsliderval');
            Q = round(src.Value);
            stxt.String = num2str(Q);     %update slider text
            %figure axes and update plot
            figax = findobj(src.Parent,'Tag','PlotFigAxes'); 
            cstPlot(obj,figax,Q)
        end
%%
        function setYaxisLimits(obj,ax)
            %set the Y axis limits so they do not change when plot updated
            dst = obj.CSTmodel;
            if ~isprop(dst,'MeanTideLevel')
                dst = activatedynamicprops(dst);
            end
            z = dst.MeanTideLevel;
            a = dst.TidalVelAmp;
            d = dst.HydDepth;
            U = dst.TidalVelAmp;
            v = dst.RiverVel;
            
            lim1 = floor(min(z-d,[],'All'));
            lim2 = ceil(max(z+a,[],'All'));
            lim3 = floor(min(v,[],'All'));
            lim4 = ceil(max(U,[],'All'));
            yyaxis left                          %fix y-axis limits
            ax.YLim = [lim1,lim2];
            yyaxis right
            ax.YLim = [lim3,lim4];
        end        
%%
        function [inp,rnp] = getModelParameters(obj,mobj)
            %extract the additional parameters needed to run the CSTmodel
            %from the Estuary and WaterLevels classes
            inp = getPropertiesStruct(obj);
            est = getClassObj(mobj,'Inputs','Estuary');
            wlv = getClassObj(mobj,'Inputs','WaterLevels');
            
            %inp parameters requires the following
            inp.WidthELength = est.WidthELength;    %width convergence length (m) =0 import from file
            inp.AreaELength = est.AreaELength;      %area convergence length (m)  =0 import from file
            inp.MTLatMouth = wlv.MSL0;              %mean tide level at mouth (mOD)
            inp.TidalAmplitude = wlv.TidalAmp;      %tidal amplitude (m)
            inp.TidalPeriod = wlv.TidalPeriod;      %tidal period (hr)
            inp.RiverDischarge = []                 ;%river discharge (m^3/s) +ve downstream

            %rnp parameters requires the following
            rnp.TimeInt = 1;     %time increment in analytical model (hrs)
            rnp.DistInt = 1000;  %distance increment along estuary (m)
            rnp.useObs = false;  %flag to indicate whether to use observations
            rnp.isfull = false;  %flag used in EnergyFlux for detail output
        end
%%
        function dsp1 = modelDSproperties(~)
            %define a dsproperties struct and add the model metadata
            dsp1 = struct('Variables',[],'Row',[],'Dimensions',[]);
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique

            %struct entries are cell arrays and can be column or row vectors
            %static ouput (mean tide values)
            dsp1.Variables = struct(...
                'Name',{'MeanTideLevel','TidalElevAmp','TidalVelAmp',...
                'RiverVel','HydDepth'},...
                'Description',{'Mean water level',...
                'Tidal elevation amplitude',...
                'Tidal velocity amplitude',...
                'River flow velocity',...
                'Hydraulic depth'},...
                'Unit',{'m','m','m/s','m/s','m'},...
                'Label',{'Mean water level (m)',...
                'Elevation amplitude (m)',...
                'Velocity amplitude (m/s)',...
                'Velocity (m/s)','Depth (m)'},...
                'QCflag',repmat({'model'},1,5));
            dsp1.Row = struct(...
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});
            dsp1.Dimensions = struct(...
                'Name',{'X'},...
                'Description',{'Chainage'},...
                'Unit',{'m'},...
                'Label',{'Distance from mouth (m)'},...
                'Format',{'-'});
        end
    end
%%        
        %add other functions to operate on properties as required   
end