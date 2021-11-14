classdef WaterLevels < muiPropertyUI            
%
%-------class help------------------------------------------------------===
% NAME
%   WaterLevels.m
% PURPOSE
%   Class for water level properties and handling change over time
% USAGE
%   obj = WaterLevels.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {' Tidal Amplitude (m)',...
                          ' Tidal Period (hr)',...
                          ' LW/HW amplitude ratio', ...
                          ' Mean Sea Level at t=0 (mAD)', ...
                          ' Rate of sea level rise (m/year)', ...
                          ' Number of Cycles',...
                          ' Amplitude (m)',...
                          ' Period (years)', ...
                          ' Phase (years)',...
                          ' M2 tidal amplitude (m)',...
                          ' S2 tidal amplitude (m)',...
                          ' O1 tidal amplitude (m)',...
                          ' M4 tidal amplitude (m)',...
                          ' M4 tidal phase (deg)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed
        %properties in struct used to call simple_tide but not used in asmita
        ElevPhase = 0      %phase of elevation (ie k.x) (rads)
        VelocityAmp = 0    %tidal velocity amplitude (m/s)
        VelocityPhase = 0  %phase of velocity (ie k.x+phi) (rads)        
    end
    
    properties
        TidalAmp = 0        %tidal amplitude (m)
        TidalPeriod = 12.4  %tidal period (hrs)
        LWtoHWratio = 1     %ratio of LW to HW tidal amplitude
        MSL0 = 0            %Mean Sea Level at t=0
        SLRrate = 0         %rate of slr (m/year)
        NumCycles = 0       %number of cycles to include
        CycleAmp = 0        %vector of cycle amplitudes (m)
        CyclePeriod = 0     %vector of cycle periods (yrs)
        CyclePhase = 0      %vector of cycle phases (yrs)
        M2amplitude = 0     %amplitude of M2 tidal constituent (m)
        S2amplitude = 0     %amplitude of S2 tidal constituent (m)
        O1amplitude = 0     %amplitude of O1 tidal constituent (m)
        M4amplitude = 0     %amplitude of M4 tidal constituent (m)
        M4phase = 0         %phase of M4 tidal constituent (m)
    end   

    properties (Dependent)
        AngularPeriod       %vector angular frequencies of periods
        AngularPhase        %vector angular frequencies of phases
        TidalRange          %tidal range used instead of TidalAmp in code
    end

    properties (Transient)
        MeanSeaLevel        %Mean sea level based on z0+slr*t
        HWaterLevel         %High Water Level at time t 
        LWaterLevel         %Low Water Level at time t
        dHWchange           %change over a time step at high water
        dLWchange           %change over a time step at low water
        HW0                 %high water at start of time step (updates each step)
        LW0                 %low water at start of time step (updates each step)
        dtr0                %amplitude offset at start of run due to imposed cycles
        dslr                %rate of sea level rise at time t (only changes if SLRrate is not constant)
    end    

%%   
    methods (Access=protected)
        function obj = WaterLevels(mobj)  
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
            classname = 'WaterLevels';              
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = WaterLevels(mobj);           
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj,9);  
                %add any additional manipulation of the input here
                if obj.NumCycles>1
                    if length(obj.CycleAmp)~=obj.NumCycles || ...
                                length(obj.CyclePeriod)~=obj.NumCycles || ...
                                         length(obj.CyclePhase)~=obj.NumCycles
                        warndlg('Amplitude, Period and Phase must equal number of cycles')
                        uiwait(gcf);
                    end
                end                
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end  
%%
        function setWaterLevels(mobj,robj)
            %update water levels for the defined time step. Called by
            %AsmitaModel.intialiseModelParameters and AsmitatModel.InitTimeStep
            obj = getClassObj(mobj,'Inputs','WaterLevels');
            mtime = robj.Time;  %model time in seconds from model t=0
            rnpobj = getClassObj(mobj,'Inputs','RunProperties');
            startyr = rnpobj.StartYear*mobj.Constants.y2s;  %startyear in seconds form Julian 0
            obj = newWaterLevel(obj,mtime,startyr);
            setClassObj(mobj,'Inputs','WaterLevels',obj);
        end
    end
%%        
        %add other functions to operate on properties as required   
        
%%
    methods
        function freq = get.AngularPeriod(obj)
            %get dependent property covert cycle period to frequency
            freq = 2*pi./obj.CyclePeriod/muiConstants.y2s;
            freq(isinf(freq)) = 0;
        end
%%
        function freq = get.AngularPhase(obj)
            %get dependent property covert angualr period to frequency
            freq = obj.AngularPeriod.*obj.CyclePhase*muiConstants.y2s;
            freq(isinf(freq)) = 0;
        end  
%%
        function tr = get.TidalRange(obj)
            %get tidal range from tidal amplitude
            tr = obj.TidalAmp*2;
        end
%%
        function initialWL(obj,startyr)
            %initialise the water level calculation
            %startyr is the start year of the simulation in seconds
            tr = obj.TidalRange;
            trfc = obj.LWtoHWratio;
            obj.HW0 = obj.MSL0 + tr/(1+trfc);
            obj.LW0 = obj.MSL0 - tr*trfc/(1+trfc);
            %tidal range amplitude offset at start of run due to imposed cycles
            obj.dtr0 = 0;
            if obj.CycleAmp(1)>0
                amp = obj.CycleAmp;
                om1 = obj.AngularPeriod; %get fcn converts input to s^-1
                om2 = obj.AngularPhase;
                %adjustment to tidal range at t=0
                obj.dtr0 = real(sum(amp.*exp(1i*(om1*startyr+om2))));
            end
        end
%%
        function [HW,MSL,LW] = newWaterLevels(obj,mtime,startyr)
            %array of water levels for model times, mtime and start year 
            %of simulation in seconds from Julian 0
            nrec = length(mtime);
            HW = zeros(nrec,1); MSL = HW; LW = HW;
            for i=1:nrec
                obj = newWaterLevel(obj,mtime(i),startyr);
                HW(i) = obj.HWaterLevel;
                LW(i) = obj.LWaterLevel;
                MSL(i) = obj.MeanSeaLevel;
            end
        end
%%
        function obj = newWaterLevel(obj,mtime,startyr)
            %obtain updated water levels at given time step
            % mtime - model time from model t=0 in seconds
            % startyr - start year of simulation in seconds from Julian 0
            tr = obj.TidalRange;
            trfc = obj.LWtoHWratio;
            amp = obj.CycleAmp;
            om1 = obj.AngularPeriod; %get fcn converts input to s^-1
            om2 = obj.AngularPhase;
            z0 = obj.MSL0;   
            %
            if mtime==0 || isempty(obj.HW0)
                initialWL(obj,startyr);
            end
            %
            if obj.SLRrate<0 
                %user has only defined the exponential rate, ert
                dslrvec = [0,1900,0.001];
                dslrvec(1) = -obj.SLRrate;   %only ert specified
                option = 2;  %exponential rate
            elseif length(obj.SLRrate)==3
                %user has defined [ert,yr0,dslr0] see model_slr for details
                dslrvec = obj.SLRrate; 
                option = 2;  %exponential rate
            elseif isnan(obj.SLRrate)
                option = 3;  %allow user to define slr function
            else
                dslrvec =  obj.SLRrate; 
                option = 1;  %linear rate
            end   
            y2s = muiConstants.y2s;   %factor to convert years to seconds
            realtime = startyr+mtime; %time in seconds from Julian 0            
            %time inputs to sealevelrise function are in years
            realyears = realtime/y2s; %convert to time in calender years
            [slr,obj.dslr] = sealevelrise(realyears,startyr/y2s,...
                                                        dslrvec,option);  
            
            msl = z0+slr; 
            dtr = real(sum(amp.*exp(1i*(om1*realtime+om2)))); %change in tidal range at t
            obj.MeanSeaLevel = msl;    
            %new water levels 
            obj.HWaterLevel = msl-obj.dtr0 + (tr+dtr)/(1+trfc); %0.5tr if trfc=1           
            obj.LWaterLevel = msl-obj.dtr0 - (tr+dtr)*trfc/(1+trfc);
            %difference with old water levels
            obj.dHWchange = obj.HWaterLevel-obj.HW0;
            obj.dLWchange = obj.LWaterLevel-obj.LW0;  
            %update old water levels to new values
            obj.HW0 = obj.HWaterLevel;   
            obj.LW0 = obj.LWaterLevel;
        end
%%
        function tabPlot(obj,src,mobj)
            %plot water levels on the WaterLevel tab
            msgtxt = 'Run properties have not been set';
            rnpobj = getClassObj(mobj,'Inputs','RunProperties',msgtxt);
            if isempty(rnpobj), return; end
            y2s = mobj.Constants.y2s;
            dt = rnpobj.TimeStep;
            nstep = rnpobj.NumSteps;            
            mtime = (0:dt:nstep*dt)*y2s;
            [HWL,MSL,LWL] = newWaterLevels(obj,mtime,rnpobj.StartYear*y2s);
            
            %clear any existing plot
            ht = findobj(src,'Type','axes');
            delete(ht);
            ax = axes(src,'Tag','PlotFigAxes');
            
            %create plot of water level change over the model run period
            ptime = rnpobj.StartYear+mtime/y2s;
            plot(ax,ptime,MSL,'Color',mcolor('orange'),'DisplayName','Mean water level')
            hold on
            plot(ax,ptime,HWL,'Color',mcolor('dark blue'),'DisplayName','High water level')
            plot(ax,ptime,LWL,'Color',mcolor('light blue'),'DisplayName','Low water level')
            hold off
            xlabel('Time (years)')
            ylabel('Elevation (mAD)')
            title('Water levels based on current settings')
            legend
        end
    end
end