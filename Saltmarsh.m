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
                          'Edge erosion rate (m/yr)',...
                          'Include marsh flat erosion (0/1)'};
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
        EdgeErosion=0;      %lateral erosion of marsh(m/yr)
        FlatErosion=false;  %flag to allow erosion of marsh surface (0/1)
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
        
%% ------------------------------------------------------------------------
% Static functions called by AsmitaModel
%-------------------------------------------------------------------------- 
        function ok = initialiseMarshDepthConc(mobj)
            %initialise the saltmarsh concentration table to define 
            %concentrations over the marsh/flat at a range of depths           
            msgtxt = 'Saltmarsh parameters not defined';
            obj = getClassObj(mobj,'Inputs','Saltmarsh',msgtxt);
            if isempty(obj), ok = 0; return; end
            
            eleobj = getClassObj(mobj,'Inputs','Element');
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            [aws,c0] = getMarshVerticalExchange(obj,eleobj);    
            %marsh concentration options
            mco.tsn = 14.77;          %duration of spring-neap cycle (days) 
            mco.delt = 10;            %time step (secs)  *may be sensitive
            mco.dmin = 0.05;          %minimum depth used in calculation           
            obj.MarshDepthConc = concovermarsh(obj,wlvobj,c0,aws,mco); 

            setClassObj(mobj,'Inputs','Saltmarsh',obj);
            ok = 1;
        end
%%
        function Bk = BioProduction(mobj)
            % Bk is a diagonal matrix for the total sedimentation and is the product 
            % of the species biomass and the rate(kBm) - units of 1/s.
            [obj,~,eleobj,cn] = Saltmarsh.getInputData(mobj);
            eletype = getEleProp(eleobj,'EleType');
            nele = length(eletype);
            ism = find(strcmp(eletype,'Saltmarsh'));
            nsm = length(ism);
            
            %get matrix of the biomass with a row for each element and
            %a column for each species (kg/m^2) 
            Bm = Biomass(obj,eleobj);
            
            %assign rate of biomass production to saltmarsh elements
            kBm = getEleSpeciesProduct(obj,eleobj,cn);
            
            %sBk is a vector for the total sedimentation and is the product 
            %of the species biomass and the rate(kBm) - units of 1/s.
            sBk = zeros(nele,1);
            for k=1:nsm
                sBk(ism(k))  = Bm(ism(k),:)*kBm(ism(k),:)';
            end
            Bk = diag(sBk);
        end
%%    
        function Deq = EqDepth(mobj)
            %calculate the equilibrium depth over saltmarsh elements
            %get model parameters
            [obj,wlvobj,eleobj,cn] = Saltmarsh.getInputData(mobj);
            dmx = obj.MaxSpDepth;
            dslr = wlvobj.dslr/mobj.Constants.y2s; %convert to seconds            
            eletype = getEleProp(eleobj,'EleType');            
            ism = find(strcmp(eletype,'Saltmarsh'));
            Deq = zeros(length(eletype),1); 
            %
            sm = getEleProp(eleobj,'SurfaceArea');
            vm = getEleProp(eleobj,'MovingVolume');
            depth = vm./sm;          % water depth in element
            if any(sm(:)==0), depth(sm==0) = 0; end  %trap infinity
            %vertical exchange for marsh (initial inorganic values)  
            ws = getEleProp(eleobj,'VerticalExchange');
            %compute equilibrium depth based on local concentration
            idp = depth<=max(dmx);   %element depth is within max marsh depth               
            ct = obj.MarshDepthConc; %concentration over marsh as a function of depth
            cem = interp1q(ct.Depth,ct.Concentration,depth);
            %if cem(ism)==0, depth(ism) = 0; end

            nsm = length(ism);
            for jsm = 1:nsm   
                if idp(ism(jsm))     %true when within marsh species range   
                    if depth(ism(jsm))>0
                        %water depth over marsh is within marsh species range
                        wsj = bioenhancedsettling(obj,depth(ism(jsm)),ws(ism(jsm)));
                        qm = wsj.*cem(ism(jsm))./depth(ism(jsm));       %sediment delivery to marsh 
                        Deq(ism(jsm)) = morris_eqdepth(obj,cn,qm,dslr); %equilibrium depth
                        if Deq(ism(jsm))==0, Deq(ism(jsm)) = -1; end
                    else                        
                        %depth of marsh has gone to zero
                        Deq(ism(jsm)) = 0;
                    end 
                else
                    %when depth greater than maximum species depth
                    %force equilibrium to the bare flat prism based value 
                    %in ASM_model.asmitaEqFunctions
                    Deq(ism(jsm)) = -1; 
                end
            end
        end
        
%%
        function ws = BioSettlingRate(mobj)
            %modify the vertical exchange vector for the influence of biomass
            %However, to avoid enhancing the erosional exchange, only do 
            %this if the element is accreting.  This is determined by 
            %whether gma=ve/vm <= 1. For gma>1 set ws to a very small value 
            %to limit erosion (gma <1: import and gma >1: export)
            eleobj = getClassObj(mobj,'Inputs','Element');
            vm = getEleProp(eleobj,'MovingVolume');
            ve = getEleProp(eleobj,'EqVolume');
            ws = getEleProp(eleobj,'transVertExch');
            sm = getEleProp(eleobj,'SurfaceArea');
            eletype = getEleProp(eleobj,'EleType');
            obj  = getClassObj(mobj,'Inputs','Saltmarsh');
            isEro = obj.FlatErosion;
            depth = vm./sm; % water depth in element
            gma = ve./vm;   % <1 import and >1 export
            ism = find(strcmp(eletype,'Saltmarsh'));
            nsm = length(ism);
            for k=1:nsm
                if gma(ism(k))> 1 && isEro
                    %marsh eroding (export) and erosion allowed
                    ws(ism(k)) = bioenhancedsettling(obj,depth(ism(k)),ws(ism(k)));
                elseif gma(ism(k))<= 1
                    %marsh accreting (import)
                    ws(ism(k)) = bioenhancedsettling(obj,depth(ism(k)),ws(ism(k)));
                else
                    %prevent erosion if gma>1 and erosion not allowed
                    %this is applied even if no biological production provided smflg=1
                    %so erosion of marsh elements is constrained regardless of biology.
                    %if smflg=0 this routine is not called and the specified values
                    %of ws are used
                    Tsn  = 14.77*24*3600; %duration of spring-neap cycle (secs)
                    ws(ism(k)) = ws(ism(k))/Tsn; %use a very small value so that matrix
                    %calculations in asmita.m that use W or Ws remain stable
                end
            end            
        end
%%
        function Blim = limitingBiomass(mobj)
            %maximum biomass*productivity for elements with positive biomass
            %called by Estuary.Response when model may not be initialised
            [obj,~,eleobj,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj)
                Blim = 0;
                return;
            end
            bmx = obj.MaxBiomass;
            if iscolumn(bmx), bmx = bmx'; end  %force a column vector
            Bm = Biomass(obj,eleobj);
            kBm = getEleSpeciesProduct(obj,eleobj,cn);
            Blim = diag(bmx*(kBm.*(Bm>0))');
        end
%%
        function [dv_ero,ds_ero] = getEdgeErosion(mobj,robj)
            %compute change in marsh and flat volumes due to edge erosion
            obj = getClassObj(mobj,'Inputs','Saltmarsh');            
            if obj.EdgeErosion==0
                dv_ero = 0; ds_ero = 0;
                return; 
            end
            ero = obj.EdgeErosion/mobj.Constants.y2s;
            dmx = max(obj.MaxSpDepth);
            %get transient Element properties for area and volume
            eleobj = getClassObj(mobj,'Inputs','Element');
            %get model parameters
            ve = getEleProp(eleobj,'EqVolume');            
            se = getEleProp(eleobj,'EqSurfaceArea');
            Le = getEleProp(eleobj,'Length');
            n  = getEleProp(eleobj,'TransportCoeff');
            eletype = getEleProp(eleobj,'transEleType');
            elename = getEleProp(eleobj,'EleName');
            %use network definition to find adjacency
            estobj = getClassObj(mobj,'Inputs','Estuary');
            estgraph = estobj.DispersionGraph;
            ism = strcmp(eletype,'Saltmarsh');
            dt = robj.delta;           %time step in seconds
            ds_ero = -ero*ism.*Le*dt;  %reduce area by erosion area
            dv_ero = ds_ero*dmx/2;     %reduce volume using half max marsh depth

            %find adjacent flats and apply the opposite change
            MarshNodes = matches(estgraph.Edges.EndNodes(:,2),elename(ism));
            flatNames = estgraph.Edges.EndNodes(MarshNodes,1);
            ismf = matches(elename,flatNames); %flats adjacent to marsh            
            ds_ero(ismf) = -ds_ero(ism);
            %if flats are sediment volumes eroded marsh volumes are removed
            %from the flat volume
            dv_ero(ismf) = -sign(n(ism)).*dv_ero(ism);

            %update the equilibrium values
            erove = num2cell(ve+dv_ero);
            erose = num2cell(se+ds_ero);
            [eleobj.EqVolume] = erove{:};  
            [eleobj.EqSurfaceArea] = erose{:};
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
%% ------------------------------------------------------------------------
% Static functions for plots and animations to display aspects of model setup
%--------------------------------------------------------------------------
        function EqDepthBiomassPlot(mobj)
            %examine influence of biomass production rates on equilibrium depth
            %produces three graphs and displays the resultant eq.depth
            %get input parameters             
            [obj,wlvobj,eleobj,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wlvobj), return; end
            %--------------------------------------------------------------
            % Store original values of kbm so that they can be restored
            kbm0 = obj.SpeciesProduct;     %species productivity (m2/kg/yr)
            %-------------------------------------------------------------- 
            [sm,ct] = initialiseSaltmarshModel(obj,wlvobj,eleobj,cn,mobj);
            if isempty(ct), return; end
            [~,c0] = getMarshVerticalExchange(obj,eleobj);  
            %--------------------------------------------------------------
            % Calculate variation with slr
            %--------------------------------------------------------------
            nint = 100;         %number of intervals
            minslr = 0.0001;    %starting value of slr
            deq = zeros(nint,1); slr = deq; biom = deq;
            for jd = 1:nint %x axis
                slr(jd) = minslr*jd; %rate of slr in m/yr
                [dep,~] = interpdepthload(obj,cn,sm.aws,c0,slr(jd)/cn.y2s);%uses (sm,cn,aws,qm0,dslr)
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
            ptxt = struct('kbm',sm.userkbm,'dp0',sm.dp0,'dp1',sm.dp1,...
                        'Qm0',sm.Qm0,'Qm1',Qm1,'Dslr',Dslr,'bm1',sm.bm1,...
                        'minslr',minslr,'maxslr',minslr*nint);
            bioInfluencePlot(obj,cn,slr,deq,biom,ptxt);
            % Advise state of marsh if elements defined
            marshElementCheck(obj,mobj)
        end
%%
        function BiomassDistributionPlot(mobj)
            %plot the distribution of each species on the bare flat profile
            [obj,wl,~,~] = Saltmarsh.getInputData(mobj);
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
            [obj,wlvobj,eleobj,cn] = Saltmarsh.getInputData(mobj);
            if isempty(obj) || isempty(wlvobj), return; end
            %--------------------------------------------------------------
            % Store original values of kbm so that they can be restored
            kbm0 = obj.SpeciesProduct;     %species productivity (m2/kg/yr)
            %-------------------------------------------------------------- 
            [sm,ct] = initialiseSaltmarshModel(obj,wlvobj,eleobj,cn,mobj);
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
                    sumKB = sum(sm.userkbm.*(bm.*(bm>0)))/cn.y2s;  
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
        
%%
        function setSaltmarsh(mobj)
            %initialise an empty instance of Saltmarsh (used in asm_oo2mui)
            obj = Saltmarsh;
            setClassObj(mobj,'Inputs','Saltmarsh',obj);
        end        
    end
    
%% ------------------------------------------------------------------------
% Open class functions
%--------------------------------------------------------------------------
    methods
        function tabPlot(obj,src,mobj)
            %add plot to Saltmarsh tab
            hp = findobj(src,'Type','uipanel');
            if isempty(hp)
                hp = uipanel(src,'Title','Biomass distribution','FontSize',8,...
                                'BackgroundColor',[0.96,0.96,0.96],...
                                'Position',[0.55 0 0.45 0.99]);
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
    
%% ------------------------------------------------------------------------
% Private functions called by static functions BioProduction and EqDepth
%-------------------------------------------------------------------------- 
    methods (Access=private)      
        function Bm = Biomass(obj,eleobj)
            %calculate the biomass for a given depth in a saltmarsh element
            % Bm is a matrix of the biomass with a row for each element and
            % a column for each species (kg/m^2) 
            
            %get model parameters
            vm = getEleProp(eleobj,'MovingVolume');
            sm = getEleProp(eleobj,'SurfaceArea');
            eletype = getEleProp(eleobj,'EleType');
            nele = length(eletype);
            ism = find(strcmp(eletype,'Saltmarsh'));
            nsm = length(ism);
            nsp = obj.NumSpecies;
            
            depth =vm./sm;
            Bc = morris_biocoeffs(obj);
            n1  = ones(nele,1);
            dd  = [depth depth.^2 n1];
           
            Bm  = zeros(nele,nsp);
            for k=1:nsm
                Bm(ism(k),:) = Bc*dd(ism(k),:)';
                jk = find(Bm(ism(k),:) <= 0);
                if ~isempty(jk)
                    Bm(ism(k),jk)= 0;
                end
            end
        end
%%
        function kBm = getEleSpeciesProduct(obj,eleobj,cn)
            %get the species productivity for each element
            eletype = getEleProp(eleobj,'EleType');
            nele = length(eletype);
            ism = find(strcmp(eletype,'Saltmarsh'));
            nsp = obj.NumSpecies;
            kbm  = obj.SpeciesProduct;
            temp = ones(nele,1);
            kBm  = zeros(nele,nsp);
            if ~isempty(kbm)
                kbm = kbm/cn.y2s;         %convert years to seconds
                for k = 1:nsp
                    kBm(ism,k) = kbm(k)./temp(ism);
                end
            end
        end  
%%
        function [aws,c0] = getMarshVerticalExchange(~,eleobj)
            %get the mean values of the vertical exchange of the saltmarsh 
            %(excluding bio) and tidal flat and the equilibirum concentration
            %over the marsh or tidal flat (only varies if fine and coarse 
            %fractions used)
            eletype = getEleProp(eleobj,'EleType');
            ws = getEleProp(eleobj,'VerticalExchange');
            ism = strcmp(eletype,'Saltmarsh');
            ifl = strcmp(eletype,'Tidalflat');
            ich = strcmp(eletype,'Channel');
            if any(ism) 
                awm = mean(ws(ism));
            else       %no saltmarsh elements defined
                awm = mean(ws(ifl));
            end
            awf = mean(ws(ifl)); 
            aws = [awm,awf];
            %equilibrium concentration over flats - only varies when fine/coarse fractions used
            cnc  = getEleProp(eleobj,'EqConcentration'); %use cE for flat elements if present
            if any(ifl) 
                c0 = mean(cnc(ifl)); %flat value - 'mean' to handle when more than one flat
            else
                c0 = cnc(1); %if no flats all eq.concs must be the same
            end
        end   
        
%% ------------------------------------------------------------------------
% Private functions for plots and animations to display aspects of model setup
%--------------------------------------------------------------------------
        function [sm,ct] = initialiseSaltmarshModel(obj,wlvobj,eleobj,cn,mobj)
            %Set up inputs needed by MarshFlatAnimation and EqDepthBiomass             
            newWaterLevels(wlvobj,0,0);                
            sm.dslr = wlvobj.SLRrate/cn.y2s;  %rate of sea level change (m/s)
            
            %intitialise transient properties            
            Element.initialiseElements(mobj);
            Element.setEqConcentration(mobj);
            [sm.aws,c0] = getMarshVerticalExchange(obj,eleobj);    
            
            %marsh concentration options
            mco.tsn = 14.77;          %duration of spring-neap cycle (days) 
            mco.delt = 10;            %time step (secs)  *may be sensitive
            mco.dmin = 0.01;          %minimum depth used in calculation           
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
            sm.userkbm = obj.SpeciesProduct;         %values in years
            %
            sm.Qm0 = 0.00018;       %estimate of sediment load used by Morris,2006
            if mean(sm.userkbm)< 1.0e-10*cn.y2s
                sm.Qm0 = 0.0018;    %adjustment needed if kbm very low, Morris, 2007
            end
            qm0 = sm.Qm0/cn.y2s;    %initial value of qm (s^-1)
            dp0 = morris_eqdepth(obj,cn,qm0,sm.dslr);
            if dp0<=0
                dmn = obj.MinSpDepth;          %minimum depth for each species (m)
                dmx = obj.MaxSpDepth;          %maximum depth for each species (m)
                dp0 = mean((dmn+dmx)/2);  
            end
            ct = obj.MarshDepthConc; %concentration over marsh as a function of depth
            cem = interp1q(ct.Depth,ct.Concentration,dp0);            
            wsm = bioenhancedsettling(obj,dp0,sm.aws);
            qm0 = wsm*cem/dp0;
            sm.dp0 = morris_eqdepth(obj,cn,qm0,sm.dslr);
            %--------------------------------------------------------------
            % Calculate depth, dp1, sediment loading, qm1, biomass
            % coefficients Bc, total biomass at equilibrium depth, bm1
            %--------------------------------------------------------------
            [sm.dp1,sm.qm1] = interpdepthload(obj,cn,sm.aws,c0,sm.dslr);
            sm.Bc = morris_biocoeffs(obj);
            dd1 = [sm.dp1 sm.dp1.^2 1];
            bm0 = (sm.Bc*dd1');
            sm.bm1 = sum(bm0.*(bm0>0)); %total biomass at equilibrium depth (kg.m^-2)
            sm.dmx = max(obj.MaxSpDepth); 
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
            ylabel('Equilibrium production (kg.m^-2)');
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
            title(sprintf('Equilibrium conditions as a function of sea level rise\n'));
        end
%%
        function marshElementCheck(obj,mobj)
            %if elements have been defined, check whether marsh is
            %biologically productive or not and display message
            msgtxt = 'Elements have not been defined';
            eleobj = getClassObj(mobj,'Inputs','Element',msgtxt);
            if ~isempty(eleobj)
                vm = getEleProp(eleobj,'MovingVolume');
                sm = getEleProp(eleobj,'SurfaceArea');
                eletype = getEleProp(eleobj,'EleType');
                elename = getEleProp(eleobj,'EleName');                
                ism = find(strcmp(eletype,'Saltmarsh')); 
                nsm = length(ism);
                dm = vm./sm;
                dmtxt = num2str(dm(ism),'%.3f');  
                maxdmx = max(obj.MaxSpDepth);
                smtxt = repmat("m is < maximum species depth",nsm,1);
                warntxt = 'All marshes are biologically active';
                if any(dm(ism)>maxdmx)
                    idx = dm(ism)>maxdmx;
                    smtxt(idx) = "m is > maximum species depth";
                    warntxt = sprintf('Some marshes have drowned\nElement depth too large, or species range too small');
                end
                dmtxt = strcat(elename(ism),repmat(" depth (Vm/Sm) of ",nsm,1),...
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
        end
%%
        function bioDistributionPlot(obj,y,z,biomass,depth,a,hfig)
            %plot of profile and biomass distribution for each species
            if isempty(hfig)
                hfig = figure('Name','Biomass Plot','Tag','PlotFig');
            end
            
            ax1_pos = [0.165,0.11,0.65,0.79]; % position of first axes
            ax1 = axes(hfig,'Position',ax1_pos,...
                      'XAxisLocation','bottom','YAxisLocation','left');
            
            profilePlot(obj,y,z,a,ax1); %plot tidal flat profile and HW  
            ax1.YLim = [ax1.YLim(1),a+0.1];
            
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
            ax.Position = [0.16,0.18,0.65,0.75]; %make space for slider bar
            hm = setSlideControl(obj,hfig,time,z,zHW,zHW-dmx);
            
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
            title('Saltmarsh development')
            nint = length(time);
            Mframes(nint,1) = getframe(gcf);
            Mframes(1,1) = getframe(gcf);
            for i=2:nint
                zi = z(i,:); zHWi = zHW(i)*[1 1]; zmxdepi = zmxdep(i)*[1 1]; %#ok<NASGU>
                refreshdata(hp,'caller')
                hm(1).Value = time(i);
                hm(3).String = string(string(time(i)));
                drawnow;                 
                Mframes(i,1) = getframe(gcf); %NB print function allows more control of resolution 
            end
            hold(ax,'off')
            
            obj.ModelMovie = Mframes; 
        end
%%
        function hm = setSlideControl(obj,hfig,time,z,zHW,zmxdep)
            %initialise slider to set different t values     
            invar = struct('sval',[],'smin',[],'smax',[],'size', [],...
                           'callback','','userdata',[],'position',[],...
                           'stxext','','butxt','','butcback','');            
            invar.sval = time(1);     %initial value for slider 
            invar.smin = time(1);     %minimum slider value
            invar.smax = time(end);   %maximum slider value
            invar.callback = @(src,evt)updateSMplot(obj,src,evt); %callback function for slider to use
            invar.userdata = struct('time',time,'z',z,'zHW',zHW,'zmxdep',zmxdep);
            invar.stext = 'Time =  ';   %text to display with slider value, if included
            invar.butxt =  'Save';    %text for button if included
            invar.butcback = @(src,evt)saveanimation2file(obj.ModelMovie,src,evt); %callback for button
            hm = setfigslider(hfig,invar);                   
        end  
%%
        function updateSMplot(~,src,~)
            %adjust the plot to the time selected by user
            stxt = findobj(src.Parent,'Tag','figsliderval');
            T = round(src.Value);
            stxt.String = num2str(T);     %update slider value text
            sldui = findobj(src.Parent,'Tag','figslider');
            sld = sldui.UserData;
            idx = find(sld.time<=T,1,'last');
            %figure axes and update plot
            figax = findobj(src.Parent.Children,'Type','Axes'); 
            hp = figax.Children;
            zi = sld.z(idx,:);               %#ok<NASGU>
            zHWi = sld.zHW(idx)*[1 1];       %#ok<NASGU>
            zmxdepi = sld.zmxdep(idx)*[1 1]; %#ok<NASGU>
            refreshdata(hp,'caller')
            drawnow;
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
        function [obj,wl,ele,cn] = getInputData(mobj)
            %initialise saltmarsh, water levels and constants
            msgtxt = 'Saltmarsh parameters not defined';
            obj = getClassObj(mobj,'Inputs','Saltmarsh',msgtxt);
            msgtxt = 'Water level data not defined';
            wl = getClassObj(mobj,'Inputs','WaterLevels',msgtxt);
            msgtxt = 'Elements not defined';
            ele = getClassObj(mobj,'Inputs','Element',msgtxt);
            cn = getConstantStruct(mobj.Constants);
        end
    end       
end