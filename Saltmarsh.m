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
            [obj,wlvobj,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wlvobj), return; end
            %--------------------------------------------------------------
            % Store original values of kbm so that they can be restored
            kbm0 = obj.SpeciesProduct;     %species productivity (m2/kg/yr)
            %-------------------------------------------------------------- 
            [sm,ct] = initialiseSaltmarshModel(obj,wlvobj,cn,mobj);
            if isempty(ct), return; end

            %--------------------------------------------------------------
            % Calculate variation with slr
            %--------------------------------------------------------------
            nint = 100;         %number of intervals
            minslr = 0.0001;    %starting value of slr
            deq = zeros(nint,1); slr = deq; biom = deq;
            for jd = 1:nint %x axis
                slr(jd) = minslr*jd; %rate of slr in m/yr
                [dep,~] = interpdepthload(obj,cn,sm.aws,sm.qm1,slr(jd)/cn.y2s);%uses (sm,cn,aws,qm0,dslr)
                deq(jd) = dep;
                if dep>0                
                    dd  = [dep dep.^2 1];
                    bm = sm.Bc*dd';
                    biom(jd) = sum(bm.*(bm>0));
                else
                    biom(jd)=0;
                end
            end
            %--------------------------------------------------------------
            % Restore original values of kbm
            obj.SpeciesProduct = kbm0;
            %--------------------------------------------------------------            
            % Plot results
            Dslr = sm.dslr*cn.y2s*1000;     %units of mm/year
            Qm1 = sm.qm1*cn.y2s;            %units of yr^-1            
            ptxt = struct('kbm',sm.kbm*cn.y2s,'dp0',sm.dp0,'dp1',sm.dp1,...
                        'Qm0',sm.Qm0,'Qm1',Qm1,'Dslr',Dslr,'bm1',sm.bm1,...
                        'minslr',minslr,'maxslr',minslr*nint);
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
            bioDistributionPlot(obj,y,z,biomass,depth,a,[]);
        end
%%
        function MarshFlatAnimation(mobj)
            %animation of the development of marsh from initial bare flat
            [obj,wlvobj,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wlvobj), return; end
            %--------------------------------------------------------------
            % Store original values of kbm so that they can be restored
            kbm0 = obj.SpeciesProduct;     %species productivity (m2/kg/yr)
            %-------------------------------------------------------------- 
            [sm,ct] = initialiseSaltmarshModel(obj,wlvobj,cn,mobj);
            if isempty(ct), return; end
            
            %prompt for run parameters
            answer = inputdlg({'MTL to HWL width:','No of years simulation',...
                                 'Start year','Include decay, 1=true'},...
                                 'Saltmarsh width',1,{'500','100','1900','0'});
            if isempty(answer), return; end
            width = str2double(answer{1});    
            nyears = str2double(answer{2});
            styear = str2double(answer{3})*cn.y2s;
            isdecay = logical(str2double(answer{4}));
            
            %get initial mud flat profile
            a = wlvobj.TidalAmp;
            [y,z0] = getFlatProfile(obj,a,width,100); %nint=100
            ymx = interp1(z0,y,(a-sm.dmx));
            
            %initialise run time parameters and water levels
            mtime = 0:1:nyears; 
            nint = length(mtime);
            mtime = mtime*cn.y2s;
            dt = 1*cn.y2s;
            [zHW,msl] = newWaterLevels(wlvobj,mtime,styear); 
            
            %compute saltmarsh elevations
            z = repmat(z0,nint,1);            
            hw = waitbar(0,'Running model');
            for i=2:nint  
                idep = find(z(i-1,:)<(zHW(i)-sm.dmx),1,'last');
                depth = zHW(i)-z(i-1,:);
                cz = interp1(ct.Depth,ct.Concentration,depth);
                %assume lower flat keeps pace with change in msl
                z(i,1:idep) = z(i-1,1:idep)+(msl(i)-msl(i-1)); %change tidalflat
                for j=idep+1:length(z)                    
                    bm = sm.Bc*[depth(j);depth(j)^2;1];
                    sumKB = sum(sm.kbm.*(bm.*(bm>0)));  
                    wsb = bioenhancedsettling(obj,depth(j),sm.aws);
                    if isdecay
                        %apply a linear decay in concentration across 
                        %the upperflat width (MTL to HWL)
                        yi = y(i)-ymx;
                        cz(j) = cz(j)*((width-yi)/(width-ymx));
                    end
                    %see eqn (4) and (9) inTownend et al, COE 2016 paper
                    % ie qm*D = wsb*1/T*integral(c*dt) == wsb*cz
                    dz = (wsb*cz(j)+sumKB*depth(j))*dt; %Krone's change in depth
                    dz(isnan(dz)) = 0;
                    z(i,j) = z(i-1,j)+dz;               %change to marsh                      
                end
                waitbar(i/nint)
            end
            close(hw)
            %--------------------------------------------------------------
            % Restore original values of kbm
            obj.SpeciesProduct = kbm0;
            %--------------------------------------------------------------
            time = (styear+mtime)/cn.y2s;            
            marshAnimationFigure(obj,y,z0,z,time,zHW,sm.dmx)
        end
    end
%%
    methods
        function tabPlot(obj,src,mobj)
            %add plot to Saltmarsh tab
            hp = findobj(src,'Type','uipanel');
            if isempty(hp)
                hp = uipanel(src,'Title','Biomass distribution','FontSize',8,...
                                'BackgroundColor',[0.96,0.96,0.96],...
                                'Position',[0.55 0 0.44 0.99]);
            end
            ht = findobj(hp,'Type','axes');
            delete(ht);
            
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            if isempty(wlvobj)
                a = 1;
            else
                a = wlvobj.TidalAmp;
            end

            [y,z] = getFlatProfile(obj,a,500,100); %nint=100
            Bc = morris_biocoeffs(obj);
            
            dmx = max(obj.MaxSpDepth);
            depth = 0:0.01:dmx;
            biomass = zeros(obj.NumSpecies,length(depth));
            for i=1:length(depth)
                    bm = Bc*[depth(i);depth(i)^2;1];
                    biomass(:,i) = bm.*(bm>0);
            end
            if isempty(biomass)
                warndlg('Some Saltmarsh properties have not been defined')
                return;
            end
            bioDistributionPlot(obj,y,z,biomass,depth,a,hp);
        end
    end
%%
    methods (Access=private)
        function [sm,ct] = initialiseSaltmarshModel(obj,wlvobj,cn,mobj)
            %Set up inputs needed by MarshFlatAnimation and EqDepthBiomass             
            newWaterLevels(wlvobj,0,0);                
            sm.dslr = wlvobj.SLRrate/cn.y2s;  %rate of sea level change (m/s)
            
            %intitialise transient properties            
            Element.initialiseElements(mobj);
            Element.setEqConcentration(mobj);
            [sm.aws,c0] = getVerticalExchange(obj,mobj);    
            
            %marsh concentration options
            mco.tsn = 14.77;          %duration of spring-neap cycle (days) 
            mco.delt = 10;            %time step (secs)  *may be sensitive
            mco.dmin = 0.05;          %minimum depth used in calculation           
            ct = concovermarsh(obj,wlvobj,c0,sm.aws,mco);         
            if all(ct.Concentration==0)
                ct = [];
                warndlg('Zero concentrations. Check Saltmarsh and Tidal Constituents are defined')
                return;
            end
            obj.MarshDepthConc = ct;
             %--------------------------------------------------------------
            % Get user defined value of kbm
            %--------------------------------------------------------------
            kbm0 = obj.SpeciesProduct;     %species productivity (m2/kg/yr)
            prompt = {'Enter biomass production rate (m^2kg^-1yr^-1)'};
            dlg_title = 'Input for biomass production rate';
            def = {num2str(kbm0)};
            dlg_ans = inputdlg(prompt,dlg_title,1,def);
            if isempty(dlg_ans), ct = []; return; end  
            obj.SpeciesProduct = str2num(dlg_ans{1}); %#ok<ST2NM>
            sm.kbm = obj.SpeciesProduct'/cn.y2s;         %values in seconds
            %
            sm.Qm0 = 0.00018;       %estimate of sediment load used by Morris,2006
            if mean(sm.kbm)< 1.0e-10
                sm.Qm0 = 0.0018;    %adjustment needed if kbm very low, Morris, 2007
            end
            qm0 = sm.Qm0/cn.y2s;       %initial value of qm (s^-1)
            sm.dp0 = morris_eqdepth(obj,cn,qm0,sm.dslr);
            %--------------------------------------------------------------
            % Calculate depth, dp1, sediment loading qm, biomass
            % coefficients Bc, total biomass at equilibrium depth and the
            % adjusted settling rate.
            %--------------------------------------------------------------
            [sm.dp1,sm.qm1] = interpdepthload(obj,cn,sm.aws,qm0,sm.dslr);
            sm.Bc = morris_biocoeffs(obj);
            dd1 = [sm.dp1 sm.dp1.^2 1];
            bm0 = (sm.Bc*dd1');
            sm.bm1 = sum(bm0.*(bm0>0)); %total biomass at equilibrium depth (kg.m^-2)
            sm.dmx = max(obj.MaxSpDepth); 
        end
%%
        function [aws,c0] = getVerticalExchange(~,mobj)
            %get the vertical exchange of the bare saltmarsh and tidal flat
            %and the concenctration over the marsh
            eleobj  = getClassObj(mobj,'Inputs','Element');
            eletype = getEleProp(eleobj,'EleType');
            ws = getEleProp(eleobj,'VerticalExchange');
            ism = strcmp(eletype,'Saltmarsh');
            ifl = strcmp(eletype,'Tidalflat');
            if any(ism) 
                awm = mean(ws(ism));
            else       %no saltmarsh elements defined
                awm = mean(ws(ifl));
            end
            awf = mean(ws(ifl)); 
            aws = [awm,awf];
            %concentration over marsh
            cnc  = getEleProp(eleobj,'EqConcentration');
            if any(ism) 
                c0 = mean(cnc(ism));
            else       %no saltmarsh elements defined
                c0 = mean(cnc(ifl));
            end
        end     
%%
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
        function bioDistributionPlot(obj,y,z,biomass,depth,a,hfig)
            %plot of profile and biomass distribution for each species
            if isempty(hfig)
                hfig = figure('Name','Biomass Plot','Tag','PlotFig');
            end
            
            ax1_pos = [0.165,0.11,0.67,0.79]; % position of first axes
            ax1 = axes(hfig,'Position',ax1_pos,...
                      'XAxisLocation','bottom','YAxisLocation','left');
            
            profilePlot(obj,y,z,a,ax1); %plot tidal flat profile and HW            
            
            style = {'-','-.','--',':'};
            green = mcolor('green');
            ax2 = axes(hfig,'Position',ax1_pos,'XAxisLocation','top',...
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
            hfig = figure('Name','Marsh animation','Tag','PlotFig');
            ax = axes(hfig,'Tag','PlotFigAxes');
            ax.Position = [0.16,0.16,0.65,0.75]; %make space for slider bar
            setSlideControl(obj,hfig,time,z,zHW,zHW-dmx);
            
            zmxdep = zHW-dmx;
            zi = z(1,:);   zHWi = zHW(1);   zmxdepi = zmxdep(1);
            profilePlot(obj,y,z0,zHWi,ax); %plot tidal flat profile and HW  
            hold(ax,'on')
            plot(ax,y,zi,'Color',mcolor('scarlet'),...
                             'LineStyle','-','DisplayName','Marsh profile')
            plot(ax,ax.XLim,zmxdepi*[1 1],'Color',mcolor('green'),...
                                         'DisplayName','Max species depth') 
            hp = ax.Children;                        
            hp(2).YDataSource = 'zi';
            hp(3).YDataSource = 'zHWi'; 
            hp(1).YDataSource = 'zmxdepi'; 
            ax.YLimMode = 'manual';   
            ax.YLim(2) = max([max(z,[],'all'),max(zHW)])+0.1;
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
            hold(ax,'off')
            
            obj.ModelMovie = Mframes; 
            answer = questdlg('Save animation','Saltmarsh animation','Yes','No','No');
            if strcmp(answer,'Yes')                 
                saveAnimation(obj)
            end
        end
%%
%%
        function hm = setSlideControl(obj,hfig,time,z,zHW,zmxdep)
            %intialise slider to set different t values          
            hm(1) = uicontrol('Parent',hfig,...
                    'Style','slider','Value',time(1),... 
                    'Min',time(1),'Max',time(end),'sliderstep',[1 1]/time(end),...
                    'Callback', @(src,evt)updateSMplot(obj,src,evt),...
                    'HorizontalAlignment', 'center',...
                    'Units','normalized', 'Position', [0.2,0.01,0.6,0.04],...
                    'UserData',struct('time',time,'z',z,'zHW',zHW,'zmxdep',zmxdep),...
                    'Tag','stepMovie');                
            hm(2) = uicontrol('Parent',hfig,...
                    'Style','text','String',num2str(time(1)),'Units','normalized',... 
                    'Position',[0.84,0.01,0.1,0.03],'Tag','SManimation');
            uicontrol('Parent',hfig,...
                    'Style','text','String','T = ','Units','normalized',... 
                    'Position',[0.81,0.01,0.05,0.03],'Tag','Ttxt');    
        end  
%%
        function updateSMplot(~,src,~)
            %adjust the plot to the time selected by user
            stxt = findobj(src.Parent,'Tag','SManimation');
            T = round(src.Value);
            stxt.String = num2str(T);     %update slider text
            sldui = findobj(src.Parent,'Tag','stepMovie');
            sld = sldui.UserData;
            idx = find(sld.time<=T,1,'last');
            %figure axes and update plot
            figax = findobj(src.Parent.Children,'Type','Axes'); 
            hp = figax.Children;
            zi = sld.z(idx,:);               %#ok<NASGU>
            zHWi = sld.zHW(idx)*[1 1];       %#ok<NASGU>
            zmxdepi = sld.zmxdep(idx)*[1 1]; %#ok<NASGU>
            refreshdata(hp,'caller')
            title(sprintf('Saltmarsh growth \nTime = %s years',string(sld.time(idx))))
            drawnow;
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
        function profilePlot(~,y,z,a,ax)
            %plot base tidal flat profile
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