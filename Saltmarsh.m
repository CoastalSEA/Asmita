classdef Saltmarsh < muiPropertyUI        
%
%-------class help------------------------------------------------------===
% NAME
%   Saltmarsh.m
% PURPOSE
%   Class to define saltmarsh and handle the influence on morphology
% USAGE
%   obj = Saltmarsh.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Number of saltmarsh species',...
                          'Minimum depth (m)','Maximum depth (m)',...
                          'Maximum biomass (kg/m2)',...
                          'Species productivity (m2/kg/yr)',...
                          'Settling cofficient, alpha (m/s)',...
                          'Settling exponent, beta (-)',...
                          'Minimum edge erosion rate (m/yr)',...
                          'Maximum edge erosion rate (m/yr)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        NumSpecies          %number of saltmarsh species
        MinSpDepth          %minimum depth for each species (m)
        MaxSpDepth          %maximum depth for each species (m)
        MaxBiomass          %maximum biomass for each species (kg/m2)
        SpeciesProduct      %species productivity (m2/kg/yr)
        SettlingAlpha       %coefficient for biomass enhanced settling rate (m/s)
        SettlingBeta        %exponent for biomass enhanced settling offset (-)
        MinEdgeErosion=0;   %minimum lateral erosion of marsh(m/yr) (only used with variable area)
        MaxEdgeErosion=0;   %maximum lateral erosion of marsh(m/yr) (ditto) 
    end    
    
    properties (Transient)
        MarshDepthConc      %lookup table of concentrations over marsh [depth,conc,submergence]
        ModelMovie          %most recent run of saltmarsh animation function
    end

%%   
    methods (Access=protected)
        function obj = Saltmarsh(mobj)  
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in the class properties defintion, above, or 
            %specify the PropertyType as a cell array here in the class 
            %constructor, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values 
            classname = 'Saltmarsh'; 
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = Saltmarsh(mobj);             
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
%%
        function EqDepthBiomass(mobj)
            %examine influence of biomass production rates on equilibirum depth
            %produces three graphs and displays the resultant eq.depth
            %get input parameters             
            [obj,wl,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wl), return; end

            newWaterLevels(wl,0,0);                
            dslr= wl.SLRrate/cn.y2s;  %rate of sea level change (m/s)
%--------------------------------------------------------------------------         
%TEMP ASSIGNMENT TO TEST CODE           
d50 = 0.0002;                      %median grain size (m)
c0 = 1.2e-4;                          %background concentration adjacent to marsh (-)
rhoc = c0*cn.rhos;
aws  = settling_velocity(d50,cn.g,cn.rhow,cn.rhos,cn.visc,rhoc);             
%             awm = mean(ws(ism));
%             awf = mean(ws(ifl));

aws = [0.003,0.0115];  %MARSH FLAT VERTICAL EXCHANGE RATE DIFFERENCE RESULTS IN DISCONTINUITY IN CONCENTRATION AT THE MARSH EDGE
                       %MAY BE BETTER TO USE SETTLING VELOCITY????
%--------------------------------------------------------------------------                        
            
            %intitialise transient properties
%             obj.h_ele = Element.initialiseElements(obj);
%             obj.h_ele = Element.setEqConcentration(obj);
%             vm = Element.getEleProp(mobj,'MovingVolume');
%             sm = Element.getEleProp(mobj,'MovingSurfaceArea');

            %marsh concentration options
            mco.tsn = 14.77;          %duration of spring-neap cycle (days) 
            mco.delt = 10;            %time step (secs)  *may be sensitive
            mco.dmin = 0.05;          %minimum depth used in calculation           
            obj.MarshDepthConc = concovermarsh(obj,wl,c0,aws,mco);            
            %--------------------------------------------------------------
            % Get user defined value of kbm
            %--------------------------------------------------------------
            kbm0 = obj.SpeciesProduct;     %species productivity (m2/kg/yr)
            prompt = {'Enter biomass production rate (m^2kg^-1yr^-1)'};
            dlg_title = 'Input for biomass production rate';
            def = {num2str(kbm0)};
            dlg_ans = inputdlg(prompt,dlg_title,1,def);
            if isempty(dlg_ans), return; end
            obj.SpeciesProduct = str2num(dlg_ans{1}); %#ok<ST2NM>
            kbm = obj.SpeciesProduct'/cn.y2s;         %values in seconds
            %
            Qm0 = 0.00018;       %estimate of sediment load used by Morris,2006
            if mean(kbm)< 1.0e-10
                Qm0 = 0.0018;    %adjustment needed if kbm very low, Morris, 2007
            end
            qm0 = Qm0/cn.y2s;       %initial value of qm (s^-1)
            dp0 = morris_eqdepth(obj,cn,qm0,dslr);
            %--------------------------------------------------------------
            % Calculate depth, dp1 and sediment loading qm
            %--------------------------------------------------------------
            [dp1,qm1] = interpdepthload(obj,cn,aws,qm0,dslr);
            Bc = morris_biocoeffs(obj);
            dd1 = [dp1 dp1.^2 1];
            bm0 = (Bc*dd1');
            bm1 = sum(bm0.*(bm0>0)); %total biomass at equilibrium depth (kg.m^-2)
            %--------------------------------------------------------------
            % Calculate variation with slr
            %--------------------------------------------------------------
            nint = 100;         %number of intervals
            minslr = 0.0001;    %starting value of slr
            deq = zeros(nint,1); slr = deq; biom = deq;
            for jd = 1:nint %x axis
                slr(jd) = minslr*jd; %rate of slr in m/yr
                [dep,~] = interpdepthload(obj,cn,aws,qm1,slr(jd)/cn.y2s);%uses (sm,cn,aws,qm0,dslr)
                deq(jd) = dep;
                if dep>0                
                    dd  = [dep dep.^2 1];
                    bm = Bc*dd';
                    biom(jd) = sum(bm.*(bm>0));
                else
                    biom(jd)=0;
                end
            end
            %--------------------------------------------------------------
            % Restore original vallues of kbm
            obj.SpeciesProduct = kbm0;
            %--------------------------------------------------------------            
            % Plot results
            Dslr = dslr*cn.y2s*1000;     %units of mm/year
            Qm1 = qm1*cn.y2s;            %units of yr^-1            
            ptxt = struct('kbm',kbm*cn.y2s,'dp0',dp0,'dp1',dp1,'Qm0',Qm0,'Qm1',Qm1,...
               'Dslr',Dslr,'bm1',bm1,'minslr',minslr,'maxslr',minslr*nint);
            bioInfluencePlot(obj,cn,slr,deq,biom,ptxt);
            % Advise state of marsh if elements defined
            marshElementCheck(obj,mobj)
        end
%%
        function BiomassDistributionPlot(mobj)
            %plot the distribution of each species on the bare flat profile
            [obj,wl,~] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wl), return; end

            width = inputdlg('MTL to HW width:','Saltmarsh width',1,{'500'});
            if isempty(width), return; end
            width = str2double(width);  
            a = wl.TidalAmp;
            [y,z] = getFlatProfile(obj,a,width,100); %nint=100
            Bc = morris_biocoeffs(obj);
            
            dmx = max(obj.MaxSpDepth);
            depth = 0:0.01:dmx;
            biomass = zeros(obj.NumSpecies,length(depth));
            for i=1:length(depth)
                    bm = Bc*[depth(i);depth(i)^2;1];
                    biomass(:,i) = bm.*(bm>0);
            end
            bioDistributionPlot(obj,y,z,biomass,depth,a);
        end
%%
        function MarshFlatAnimation(mobj)
            %animation of the development of marsh from initial bare flat
            [obj,wl,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wl), return; end

%--------------------------------------------------------------------------         
%TEMP ASSIGNMENT TO TEST CODE            
d50 = 0.0002;              %median grain size (m)
c0 = 1.2e-4;               %background concentration adjacent to marsh (-)
rhoc = c0*cn.rhos;
aws  = settling_velocity(d50,cn.g,cn.rhow,cn.rhos,cn.visc,rhoc);  
%--------------------------------------------------------------------------                     
            %marsh concentration options
            mco.tsn = 14.77;          %duration of spring-neap cycle (days) 
            mco.delt = 10;            %time step (secs)  *may be sensitive
            mco.dmin = 0.05;          %minimum depth used in calculation           
            ct = concovermarsh(obj,wl,c0,aws,mco);  
            
            dmx = max(obj.MaxSpDepth);                
            answer = inputdlg({'MTL to HW width:','No of years simulation',...
                                 'Start year','Include decay'},...
                                 'Saltmarsh width',1,{'500','100','0','0'});
            if isempty(answer), return; end
            width = str2double(answer{1});    
            nyears = str2double(answer{2});
            styear = str2double(answer{3})*cn.y2s;
            isdecay = logical(str2double(answer{4}));
            a = wl.TidalAmp;
            [y,z0] = getFlatProfile(obj,a,width,100); %nint=100
            ymx = interp1(z0,y,(a-dmx));
            Bc = morris_biocoeffs(obj);
            kbm = obj.SpeciesProduct'/cn.y2s;
            
            mtime = 0:1:nyears; 
            nint = length(mtime);
%             if rem(nint,2)~=0
%                 nint=nint+1;
%                 mtime(nint) = nyears+1;
%             end
            mtime = mtime*cn.y2s;
            dt = 1*cn.y2s;

            zHW = newWaterLevels(wl,mtime,styear); 
            z = repmat(z0,nint,1);            
            hw = waitbar(0,'Running model');
            for i=2:nint  
                idep = find(z(i-1,:)<(zHW(i)-dmx),1,'last');
                z(i,1:idep) = z(i-1,1:idep)+(zHW(i)-zHW(i-1));
                for j=idep+1:length(z)                    
                    depth = zHW(i)-z(i-1,:);
                    bm = Bc*[depth(j);depth(j)^2;1];
                    sumKB = sum(kbm.*(bm.*(bm>0)));  
                    wsb = bioenhancedsettling(obj,depth(j),[aws,aws]);
                    cz = interp1(ct.Depth,ct.Concentration,depth(j));
                    if isdecay
                        yi = y(i)-ymx;
                        cz = cz*((width-yi)/(width-ymx));
                    end
                    dz = (wsb*cz+sumKB)*depth(j)*dt;
                    dz(isnan(dz)) = 0;
                    z(i,j) = z(i-1,j)+dz;
                end
                waitbar(i/nint)
            end
            close(hw)
            time = (styear+mtime)/cn.y2s;            
            marshAnimationFigure(obj,y,z0,z,time,zHW,dmx)
        end
    end
%%
    methods (Access=private)
        function bioInfluencePlot(obj,cn,slr,deq,biom,tx)
            %plot marsh concentrations, submergence, response to SLR
            dpm = obj.MarshDepthConc.Depth;
            cem = obj.MarshDepthConc.Concentration;
            smt = obj.MarshDepthConc.Submergence;
            
            dmx = max(obj.MaxSpDepth)+0.1;
            idxmarsh = find(dpm>=dmx,1,'first');
            cemm = cem(1:idxmarsh);
            smtm = smt(1:idxmarsh);
            dpmm = dpm(1:idxmarsh);
            
            figure('Name','Biomass Plot','Tag','PlotFig');
            %--------------------------------------------------------------
            subplot(2,2,1)                     %flat and marsh - top left
            yyaxis left
            plot(dpm,cem*cn.rhos);
            xlabel('Depth over marsh/flat surface (m)');
            ylabel('Concentration over flat (kgm^-3)');
            ylim1 = get(gca,'Ylim');
            yticks(0:ylim1(2)/10:ylim1(2));
            
            yyaxis right
            plot(dpm,smt,'-.');
            ylabel('Proportion of time submerged');
            ylim([0,1])
            yticks(0:0.1:1);            
            title(sprintf('Average concentrations over intertidal flat\n'));
            
            %--------------------------------------------------------------
            subplot(2,2,3)                     %marsh only - bottom left
            yyaxis left
            plot(dpmm,cemm*cn.rhos);
            xlim([0,dmx])
            xlabel('Depth over marsh surface (m)');
            ylabel('Concentration over marsh (kgm^-3)');
            ylim1 = get(gca,'Ylim');
            yticks(0:ylim1(2)/10:ylim1(2));
            
            yyaxis right                 
            plot(dpmm,smtm,'-.');
            ylabel('Proportion of time submerged');
            ylim2 = get(gca,'Ylim');
            yticks(0:ylim2(2)/10:ylim2(2));            
            title(sprintf('Average concentrations over marsh\n'));
            
            %--------------------------------------------------------------
            subplot(2,2,[2,4])                 %slr plot - right side
            yyaxis left
            plot(slr,deq,'b-');
            ax = gca;
            ax.XScale = 'log';
            xlabel('Rate of sea level rise (m/yr)');
            ylabel('Equilibrium depth (m)'); 
            ylim([0,1]);
            ylim1 = get(gca,'Ylim');
            yticks(0:0.2:ylim1(2));
            line('XData',[tx.Dslr/1000 tx.Dslr/1000],'YData',[0 tx.dp1],'Color','r','Linestyle','--');
            line('XData',[tx.minslr tx.Dslr/1000],'YData',[tx.dp1 tx.dp1],'Color','r','Linestyle','--');           
            
            yyaxis right
            plot(slr,biom,'-.');
            ylabel('Equilibirum production (kg.m^-2)');
            ylim2 = get(gca,'Ylim');
            yticks(0:0.2:ylim2(2));
            line('XData',[tx.Dslr/1000 tx.maxslr],'YData',[tx.bm1 tx.bm1],'Color','r','Linestyle',':');
            line('XData',[tx.Dslr/1000 tx.Dslr/1000],'YData',[0 tx.bm1],'Color','r','Linestyle',':');
            
            out_str1 = sprintf('Default sediment loading, qm = %5.3d yr-1 and slr = %g mm/yr give an equilibrium depth of %5.3f m',...
                tx.Qm0,tx.Dslr,tx.dp0);
            out_str2 = sprintf('Concentrations over marsh give a sediment loading, qm = %5.3d yr-1, and a depth of %5.3f m',...
                tx.Qm1,tx.dp1);
            out_str3 = sprintf('Species productivity values used: %g %g %g',tx.kbm);
            out_str = sprintf('%s\n%s\n%s',out_str1,out_str2,out_str3);
            annotation('textbox','String',out_str,'FitBoxToText','on',...
                'HorizontalAlignment','left', 'VerticalAlignment','middle',...
                'Position',[0.6,0.8,0.28,0.06],'LineStyle','none');            
            title(sprintf('Equilibrium conditons as a function of sea level rise\n'));
        end
%%
        function marshElementCheck(obj,mobj)
            %if elements have been defined, check whether marsh is
            %biologically productive or not and display message
            msgtxt = 'Elements have not been defined';
            eleobj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if ~isempty(eleobj)
                vm = getEleProp(eleobj,'MovingVolume');
                sm = getEleProp(eleobj,'MovingSurfaceArea');
                eletype = getEleProp(eleobj,'EleType');
                elename = getEleProp(eleobj,'EleName');                
                ism = find(strcmp(eletype,'Saltmarsh')); 
                nsm = length(ism);
                dm = vm./sm;
                dmtxt = string(dm(ism));  
                maxdmx = max(obj.MaxSpDepth);
                smtxt = repmat("m is < maximum species depth",nsm,1);
                warntxt = 'All marshes are biologically active';
                if any(dm(ism)>maxdmx)
                    idx = dm(ism)>maxdmx;
                    smtxt(idx) = "m is > maximum species depth";
                    warntxt = sprintf('Some marshes have drowned\nElement depth too large, or species range too small');
                end
                dmtxt = strcat(elename(ism),repmat(" depth of ",nsm,1),...
                    dmtxt,smtxt);
                msgtxt = [sprintf('Maximum species depth is %.3g;',maxdmx);...
                                                       dmtxt;warntxt];
                msgbox(msgtxt,'Saltmarsh depths','help')
            end           
        end        
        
%%
        function [y,z] = getFlatProfile(~,a,width,nint)
            %bare flat profile based on Friedrichs tidal equilibrium form
            Ls = width/(pi/2);
            y = 0:width/nint:width;
            z = a*sin(y/Ls);
%             z1 = a*((y-Ls)/Ls).*(y<=Ls);
%             z2 = a+a*sin((y-Ls)/Ls).*(y>Ls);
%             z = z1+z2;
        end
%%
        function bioDistributionPlot(obj,y,z,biomass,depth,a)
            %plot of profile and biomass distribution for each species
            ax1 = profilePlot(obj,y,z,a); %plot tidal flat profile and HW            
            
            style = {'-','-.','--',':'};
            green = mcolor('green');
            
            ax1_pos = ax1.Position; % position of first axes
            ax2 = axes('Position',ax1_pos,'XAxisLocation','top',...
                'YAxisLocation','right','Color','none');
            ax2.XDir = 'reverse';
            ax2.YLim = ax1.YLim;
            ax2.YTickLabel = string(a-ax1.YTick);
            zd = a-depth;
            line(biomass(1,:),zd,'Parent',ax2,'Color',green,...
                              'LineStyle','-','DisplayName','Species 1')
            hold on            
            for j=2:obj.NumSpecies
                spectxt = sprintf('Species %d',j);
                line(biomass(j,:),zd,'Parent',ax2,'Color',green,...
                              'LineStyle',style{j},'DisplayName',spectxt)          
            end
            hold off            
            xlabel('Biomass (kg/m^2)')
            ylabel('Depth')
            
            legend(ax2,'Location','east')
        end
%%
        function marshAnimationFigure(obj,y,z0,z,time,zHW,dmx)
            %animation of marsh surface
            zmxdep = zHW-dmx;
            zi = z(1,:);   zHWi = zHW(1);   zmxdepi = zmxdep(1);
            ax1 = profilePlot(obj,y,z0,zHWi); %plot tidal flat profile and HW  
            hold(ax1,'on')
            plot(ax1,y,zi,'Color',mcolor('scarlet'),...
                             'LineStyle','-','DisplayName','Marsh profile')
            plot(ax1,ax1.XLim,zmxdepi*[1 1],'Color',mcolor('green'),...
                                         'DisplayName','Max species depth') 
            hp = ax1.Children;                        
            hp(2).YDataSource = 'zi';
            hp(3).YDataSource = 'zHWi'; 
            hp(1).YDataSource = 'zmxdepi'; 
            ax1.YLimMode = 'manual';   
            ax1.YLim(2) = max([max(z,[],'all'),max(zHW)])+0.1;
            title(sprintf('Saltmarsh growth \nTime = %s years',string(time(1))))
            nint = length(time);
            Mframes(nint,1) = getframe(gcf);
            Mframes(1,1) = getframe(gcf);
            for i=2:nint
                zi = z(i,:); zHWi = zHW(i)*[1 1]; zmxdepi = zmxdep(i)*[1 1]; %#ok<NASGU>
                refreshdata(hp,'caller')
                title(sprintf('Saltmarsh growth \nTime = %s years',string(time(i))))
                drawnow;                 
                Mframes(i,1) = getframe(gcf); %NB print function allows more control of resolution 
            end
            hold(ax1,'off')

            answer = questdlg('Save animation','Saltmarsh animation','Yes','No','No');
            if strcmp(answer,'Yes')
                obj.ModelMovie = Mframes;  
                saveAnimation(obj)
            end
        end
%%        
        function saveAnimation(obj)
            %save an animation plot created by PlotFig.newAnimation    
            if isempty(obj.ModelMovie)
                warndlg('No movie has been created. Create movie first');
                return;
            end
            [file,path] = uiputfile('*.mp4','Save file as','moviefile.mp4');
            if file==0, return; end
            v = VideoWriter([path,file],'MPEG-4');
            open(v);            
            warning('off','MATLAB:audiovideo:VideoWriter:mp4FramePadded')
            writeVideo(v,obj.ModelMovie);
            warning('on','MATLAB:audiovideo:VideoWriter:mp4FramePadded')
            close(v);
        end         
%%
        function ax = profilePlot(~,y,z,a)
            %plot base tidal flat profile
            hfig = figure('Name','Biomass Plot','Tag','PlotFig');
            ax = axes(hfig,'XAxisLocation','bottom','YAxisLocation','left');

            plot(ax,y,z,'Color','k','DisplayName','Tidal flat')
            hold on
            plot(ax,ax.XLim, a*[1 1],'Color','b','DisplayName','High water')
            hold off
            xlabel('Distance (m)')
            ylabel('Elevation (mOD)')
            legend(ax,'Location','southeast')
        end
    end
 %%
    methods (Static,Access=private)
        function [obj,wl,cn] = getInputData(mobj)
            %initialise saltmarsh, water levels and constants
            msgtxt = 'Saltmarsh parameters not yet defined';
            obj = getClassObj(mobj,'Inputs','Saltmarsh',msgtxt);
            msgtxt = 'Water level data not defined';
            wl = getClassObj(mobj,'Inputs','WaterLevels',msgtxt);
            cn = getConstantStruct(mobj.Constants);
        end
    end       
end