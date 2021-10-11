classdef ASMinterface < handle
    % Abstract class for ASM_model. Has all functions needed to run model
    % all of which can be overloaded in ASM_model
    %
    %----------------------------------------------------------------------
    % AUTHOR
    % Ian Townend
    %
    % COPYRIGHT
    % CoastalSEA, (c) 2021
    %----------------------------------------------------------------------
    %
    properties (Transient)
        eqScaling   %scaling of equilibirum relative to initial value
        export      %summation of net export over run (for mass balance)
        intervent   %summation of interventions over run (for mass balance)
        dWLvolume   %summation of volume changes due to water level change
        SedMbal     %sediment mass balance updated at each time step
        WatMbal     %water mass balance updated at each time step
        UniqueYears   %array of years for imposed changes
        AnnualChange  %array of imposed volume and area changes for each element
        DQ          %matix of dispersion and advection updated each time step
        dqIn        %vector of input dispersion and advection updated at each time step
    end

    methods (Abstract,Static)
        setASM_model(mobj)
    end
    
    methods (Static)
        obj = UserPrismCoeffs(mobj);    %external function

%%
        function asmitaVolumeChange(mobj,robj)
            %calculate the changes in volume and concentration over a timestep
            %and check the mass balance
            eleobj = getClassObj(mobj,'Inputs','Element');
            %get changes in water level at high and low water
            dwl = getEleProp(eleobj,'EleWLchange');
            nele = length(eleobj);
            %get model parameters
            vm = getEleProp(eleobj,'MovingVolume');
            vf = getEleProp(eleobj,'FixedVolume');
            ve = getEleProp(eleobj,'EqVolume');
            vb = getEleProp(eleobj,'BioProdVolume');
            sm = getEleProp(eleobj,'MovingSurfaceArea');
            n  = getEleProp(eleobj,'TransportCoeff');
            cb = getEleProp(eleobj,'BedConcentration');            
            
            %get concentration before updating volumes
            conc = ASM_model.asmitaConcentrations(mobj);

            [B,dd] = ASM_model.BddMatrices(mobj);
            Gam  = (ve./vm).^n;
            %Gam(vm==0) = 1;
            %trap divide by zero if vm=0
            % vis0 = find(vm<=0);
            % if ~isempty(vis0)
            %     Gam(vis0)=0;
            % end
            %
            % morphological change in volume
            dvf  = (B*Gam-dd)./cb*robj.delta;
            %adjust for biomass change in volume
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncSaltmarsh
                Bk = Saltmarsh.BioProduction(mobj);
                dvb = diag(sign(n))*(Bk*vm)*robj.delta;
            else
                dvb = zeros(size(dvf));
            end
            %
            % if ~isempty(vis0)
            %     dvf(vis0) = 0; %force dvf to zero if vm=0
            % end
            
            % change in water volume (only applies to water volumes)
            dvm = diag(n>0)*(sm.*dwl);
            %
            vm = vm + dvm + dvf - dvb; %total change (moving surface)
            vf = vf + dvf - dvb;       %morphological change (fixed surface)
            vb = vb + dvb;             %saltmarsh organic sedimentation
            
            %check that elements have not infilled
            if any(vm(:)<=0)      %when elements go to zero retain small
                vf(vm<=0) = 999;  %value to prevent matrix becoming poorly 
                vm(vm<=0) = 999;  %conditioned          
            end                      
              
            % assign results
            for i=1:nele
                eleobj(i).MovingVolume = vm(i);
                eleobj(i).FixedVolume = vf(i);
                eleobj(i).EleConcentration = conc(i);
                eleobj(i).BioProdVolume = vb(i);
            end
            %mass balance check
            ASM_model.MassBalance(mobj,robj);
            %diagnostic message to command window during runtime
%             t = robj.Time/mobj.Constants.y2s;
%             smb = obj.SedMbal; wmb = obj.WatMbal;
%             sprintf('t = %g; Sediment Balance %g; Water Balance %g',t,smb,wmb)
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
        
%%
        function [B,dd] = BddMatrices(mobj)
            %get model parameters
            eleobj = getClassObj(mobj,'Inputs','Element');
            sm = getEleProp(eleobj,'MovingSurfaceArea');
            ws = getEleProp(eleobj,'transVertExch');
            n = getEleProp(eleobj,'TransportCoeff');
            cE = getEleProp(eleobj,'EqConcentration');
            obj  = getClassObj(mobj,'Inputs','ASM_model');
            DQ = obj.DQ;  %values defined using setDQmatrix
            dqIn = obj.dqIn;
            %
            % Diagonal matrix for vertical exchange and surface areas
            W  = diag(ws);
            Sm = diag(sm);
            %
            % Diagonal matrix for: 
            M = diag(sign(n));   %sign to define volumes as wet or dry
            CE = diag(cE);       %eq conc by element type
            I = eye(length(sm)); %unit matrix
            %
            B    = CE*M*Sm*W*(I-(DQ+Sm*W)\Sm*W);
            dd   = CE*M*Sm*W*((DQ+Sm*W)\dqIn);
            %full equation is as follows(DQ and dqIn replace terms in brackets)
            % B   = cE*M*Sm*Wz*(I-(D+Q+Qs+Qtp+Sm*Wz)\Sm*Wz);
            % dd  = cE*M*Sm*Wz*((D+Q+Qs+Qtp+Sm*Wz)\(dExt+qtpIn+kCeI.*qIn+qsIn));          
        end
        
%% 
        function setDQmatrix(mobj,offset)
            %set the DQ, dqIn and conc properties to be used for a time step
            %set up dispersion and advection matrices using dynamic
            %ReachGraph and FlowGraph if IncDynamicElements is true
            %otherwise use Estuary-Dispersion and Advection-Flow properties
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj  = getClassObj(mobj,'Inputs','Element');
            cE = getEleProp(eleobj,'EqConcentration');
            kCeI = River.getRiverProp(mobj,'tsRiverConc')./cE;
            [D,dExt] = Estuary.getDispersion(mobj);
            [Q,qIn,~] = Advection.getAdvectionFlow(mobj,'River');
            [Qtp,qtpIn,~] = Advection.getAdvectionFlow(mobj,'Qtp');
            [Qs,qsIn,~] = Advection.getAdvectionFlow(mobj,'Drift');

            %to set eqCorV in Element.setEleAdvOffsets need to only
            %include some of the advections based on conditions set
            %and offset = RunConditions.Adv2Offset
            %whereas runtime calls use all advections included in run
            %and offset = RunConditions.Adv2Inc
            switch offset
                case 'flow+drift'    %include flow and drift
                    obj.DQ = D+Q+Qtp+Qs;
                    obj.dqIn = dExt+kCeI.*qIn+qtpIn+qsIn;
                case 'flow'          %include flow only
                    obj.DQ = D+Q+Qtp;
                    obj.dqIn = dExt+kCeI.*qIn+qtpIn;
                case 'drift'         %include drift only
                    obj.DQ = D+Qs;
                    obj.dqIn = dExt+qsIn;
                otherwise            %use no advections (offset = 'none')
                    obj.DQ = D;
                    obj.dqIn = dExt;
            end   

            setClassObj(mobj,'Inputs','ASM_model',obj);
        end

%%
        function setVertExch(mobj,robj)
            %get model parameters
            eleobj = getClassObj(mobj,'Inputs','Element');
            nele = length(eleobj);
            vm = getEleProp(eleobj,'MovingVolume');
            Vo = getEleProp(eleobj,'InitialVolume');
            ve = getEleProp(eleobj,'EqVolume');
            cb = getEleProp(eleobj,'BedConcentration');
            n  = getEleProp(eleobj,'TransportCoeff');
            ws = getEleProp(eleobj,'VerticalExchange');
            %adjust vertical exchange if saltmarsh included
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncSaltmarsh
                ws = Saltmarsh.BioSettlingRate(mobj);
            end
            %reset transVertExch to initial value (or as modified by bio)
            %because ws used to calculate B and dd to get unconstrained dvf
            for i=1:nele
                eleobj(i).transVertExch = ws(i);
            end

            %unadjusted morphological change in volume 
            [B,dd] = ASM_model.BddMatrices(mobj);
            Gam  = (ve./vm).^n;
            dvf  = (B*Gam-dd)./cb*robj.delta;
            
            %check if any erosion is greater than available sediment
            vf = getEleProp(eleobj,'FixedVolume');            
            noero = ~getEleProp(eleobj,'Erodible'); %ie if Erodible false            
            %when noero true and SedChange>0 ie erodiing
            %correction=1 if SedAvailable>SedChange
            %correction=0 if SedAvailable<=0
            %correction=SedAvailable/SedChange if SedAvailable<=SedChange
            SedAvailable = sign(n).*(Vo-vf);
            SedChange = sign(n).*dvf;
            for j=1:nele
                correction = 1;                
                if noero(j) && SedChange(j)>0
                    if SedAvailable(j)<=0 
                        correction = 0;  
                    elseif SedAvailable(j)<SedChange(j)
                        correction = SedAvailable(j)/SedChange(j);    
                    end
                end
                %update vertical exchange based on any correction needed
                eleobj(j).transVertExch = ws(j)*correction;                
            end
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
%%
        function asmitaEqFunctions(mobj)            
            %function to define the equilibrium volume for a given prism
            %overloaded in ASM_model
            eleobj = getClassObj(mobj,'Inputs','Element');
            ecpobj = getClassObj(mobj,'Inputs','EqCoeffProps');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            
            etypalpha = ecpobj.alpha;
            etypbeta = ecpobj.beta;
            eqType  = ecpobj.eqtype;
            %
            %Available eletypes are defined in mobj.GeoType
            eletype = getEleProp(eleobj,'transEleType');
            %
            prism = Reach.getReachEleProp(mobj,'UpstreamPrism');
            HWL = Reach.getReachEleProp(mobj,'HWlevel');
            LWL = Reach.getReachEleProp(mobj,'LWlevel');
            %
            %equilibirum surface area taken as intial area (ie fixed)
            EqSA = getEleProp(eleobj,'InitialSurfaceArea');
            %Equilibirum depth over marsh elements
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncSaltmarsh
                Deq = Saltmarsh.EqDepth(mobj); 
            else
                Deq = zeros(size(eletype));
            end
            %assign equilibrium volume to each element
            for i=1:length(eleobj)
                alpha = etypalpha.(eletype{i});
                beta = etypbeta.(eletype{i});
                isTReq = logical(eqType.(eletype{i}));
                eleobj(i).EqSurfaceArea = EqSA(i);
                switch eletype{i}
                    case 'Saltmarsh'
                        if Deq(i)>0
                            eleobj(i).EqVolume = EqSA(i)*Deq(i);
                        else
                            eleobj(i).EqVolume = alpha*prism(i)^beta;
                        end
                    otherwise
                        if ~isTReq %appplies to any element type (eg tidalflat)
                            eleobj(i).EqVolume = alpha*(HWL(i)-LWL(i))^beta;
                        else
                            eleobj(i).EqVolume = alpha*prism(i)^beta;
                        end
                end
            end
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
        
%%
        function asmitaEqScalingCoeffs(mobj)
            %function to define scaling for the equilibrium relative to
            %the initial condition          
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj = getClassObj(mobj,'Inputs','Element');
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            
            vm0 = getEleProp(eleobj,'InitialVolume');
            ASM_model.asmitaEqFunctions(mobj);
            ve = getEleProp(eleobj,'EqVolume');
            if rncobj.ScaleValues
                obj.eqScaling(:,1) = vm0./ve;
            else
                obj.eqScaling(:,1) = ones(size(vm0));
            end
            setClassObj(mobj,'Inputs','ASM_model',obj);
        end
           
%%
        function conc = asmitaConcentrations(mobj,varargin)
            %get concentration taking acount of the horizontal exchanges
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj = getClassObj(mobj,'Inputs','Element');
            
            vm = getEleProp(eleobj,'MovingVolume');
            ve = getEleProp(eleobj,'EqVolume');
            sm = getEleProp(eleobj,'MovingSurfaceArea');
            %
            ws = getEleProp(eleobj,'transVertExch');
            n = getEleProp(eleobj,'TransportCoeff');
            %
            if isempty(varargin)
                DQ = obj.DQ;  %values defined using setDQmatrix
                dqIn = obj.dqIn;
            else
                DQ = varargin{1};
                dqIn = varargin{2};
            end
            %
            cE = getEleProp(eleobj,'EqConcentration');
            %
            % Diagonal matrix for vertical exchange and surface areas
            W  = diag(ws);
            Sm = diag(sm);
            Gam  = (ve./vm).^n;
            %
            % actual concentration taking account of flow field.
            conc = (DQ+Sm*W)\(Sm*W*Gam+(dqIn)).*cE;
            %once DQ and dqIn are expanded, this is equivalent to:
            % cnm = (D+Q+Qs+Qtp+Sm*Wz)\(Sm*Wz*Gam+(dE+qtpin+kCeI.*qEin+qSin))*cE;
        end
        
 %%      
         function obj = clearMassBalance(mobj)
             obj = getClassObj(mobj,'Inputs','ASM_model');
             obj.export = 0;
             obj.intervent = 0;
             obj.dWLvolume = 0;
             obj.SedMbal = 0;
             obj.WatMbal = 0;
         end
         
%%
        function MassBalance(mobj,robj)
            obj = getClassObj(mobj,'Inputs','ASM_model');
            if robj.iStep<=1
                obj = ASMinterface.clearMassBalance(mobj);
            end
            % Load data
            delta = robj.delta;                    %time step in seconds
            time = robj.Time;                      %elapsed time in seconds

            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            eleobj = getClassObj(mobj,'Inputs','Element');
            
            eletype = getEleProp(eleobj,'transEleType');
            cE = getEleProp(eleobj,'EqConcentration');
%             cES = estobj.EqConcCoarse;
%             cEM = estobj.EqConcFine;           
            cb = getEleProp(eleobj,'BedConcentration');
            cI = River.getRiverProp(mobj,'tsRiverConc');
            conc = getEleProp(eleobj,'EleConcentration');
            n = getEleProp(eleobj,'TransportCoeff');
            M = diag(sign(n));
            wM = diag(n>0);
            
            Vo = getEleProp(eleobj,'InitialVolume');
            vf = getEleProp(eleobj,'FixedVolume');
            vm = getEleProp(eleobj,'MovingVolume');
            vb = getEleProp(eleobj,'BioProdVolume');
            dV = Interventions.getIntProp(mobj,'transVolChange');
            [~,dExt] = Estuary.getDispersion(mobj); %uses ReachGraph
            [~,qIn,qOut] = Advection.getAdvectionFlow(mobj,'River'); %water flux in m3/s
            [~,qtpIn,qtpOut] = Advection.getAdvectionFlow(mobj,'Qtp'); %water flux in m3/s
            [~,qsIn,qsOut] = Advection.getAdvectionFlow(mobj,'Drift'); %sed.flow in m3/s

            % Calculate sediment mass balance (+ve is import to water column)
            impriv = sum(qIn.*cI)*time;          %river import
            impqtp = sum(qtpIn.*cE)*time;        %tidal pumping marine import
            impsed = sum(qsIn.*cE)*time;         %littoral drift import
            import = impriv+impqtp+impsed;       %import to system
            dexp1  = sum(qOut.*conc);            %advective outflow
            dexp2  = sum(qtpOut.*conc);          %tidal pumping outflow
            dexp3  = sum(qsOut.*conc);           %littoral drift export
            dexp4  = sum(dExt.*(conc-cE));       %diffusion exchange with env (coarse)
            dexp   = (dexp1+dexp2+dexp3+dexp4)*delta;  %sum over a time step
            obj.export = obj.export+dexp;        %export to marine
            bed    = sum(M*(Vo-vf).*cb);         %bed changes
            obj.intervent = obj.intervent+sum(dV.*cb);  %user defined changes
            biosed = sum(vb.*cb);                %saltmarsh organic sedimentation
            %balance of cumulative changes from t=0 
            obj.SedMbal = import+biosed-bed-obj.export-obj.intervent;
            
            %take account of mixed sediment concentrations
%             if cES~=cEM
%                 obj.SedMbal = ASM_model.SandMudCorrection(mobj,robj);
%             end

            %calculate water mass balance (+ve is increase in water volume)             
            %change in water volume due to change in water level at high water; 
            idx = strcmp(eletype,'Tidalflat');
            if all(n(idx)<0)     %all tidal flats defined as sediment volumes
                                 %only use elements defined as water volumes   
                wetArea = wM*getEleProp(eleobj,'MovingSurfaceArea');                
            elseif any(n(idx)<0) %some tidal flats defined as sediment volumes
                wetArea = NaN;
            else                 %all tidal flats defined as water volumes                
                 wetArea = Reach.getReachProp(mobj,'HWarea');
            end                          
            dHW = wlvobj.dHWchange; 
            damp = getEleProp(eleobj,'TidalDamping');
            idr = Reach.getReachProp(mobj,'ReachChannelID');
            obj.dWLvolume = obj.dWLvolume+sum(wetArea.*(dHW*damp(idr)));
            
            %change in water volume due to volume changes (water + sediment)
            sLW = Reach.getReachProp(mobj,'LWarea');
            tr = wlvobj.HWaterLevel-wlvobj.LWaterLevel;            
            dtr = (tr-wlvobj.TidalRange)*damp;
            %wM limits the elements included to those that are water volumes
            watervols = wM*(vm-Vo);  %change in total volumes
            morphvols = wM*(Vo-vf);  %change in morphological volume
            %total change including change in prism over channel elements  
            estvol = sum(watervols)+sum(morphvols)+sum(sLW.*dtr(idr));  
            
            %check balance
            obj.WatMbal = estvol-obj.dWLvolume;

            setClassObj(mobj,'Inputs','ASM_model',obj);            
        end
        
%%
        function MassBalance = SandMudCorrection(mobj,robj)
            obj = getClassObj(mobj,'Inputs','ASM_model');
            MassBalance = obj.SedMbal;
%             delta = robj.delta;
%             idx = (strcmp(eletype,'Tidalflat'));
%             idx = idx+(strcmp(eletype,'Saltmarsh'));
%             idx = idx+(strcmp(eletype,'Storage'));
%             
%             deltamud = D(D>0).*idx;%only correct if there are no external exchanges to these element types
%             sm = getEleProp(eleobj,'MovingSurfaceArea');
%             ve = getEleProp(eleobj,'EqVolume');
%             n = getEleProp(eleobj,'TransportCoeff');
%             ws = getEleProp(eleobj,'transVertExch');
%             
%             [DQ,dqIn] = ASM_model.getDQmatrix(mobj,0);
%             W  = diag(ws);
%             Sm = diag(sm);
%             Gam  = (ve./vm).^n;
%             %
%             % actual concentration taking account of flow field.
%             cS = (DQ+Sm*W)\(Sm*W*Gam+(dqIn)).*cES;
%             cM = (DQ+Sm*W)\(Sm*W*Gam+(dqIn)).*cEM;
%             
%             disp = estobj.Dispersion;
%             
%             idr = find(disp(:,(idx>0))>0);
%             idc = find(idx>0);
%             
%             wsS = sum(ws(~idx))/sum(~idx);
%             wsM = ws.*idx;
%             
%             deltasand = deltamud.*(wsM./wsS);
%             
%             dS = disp;
%             dS(idr,idc) = deltasand(idx>0);
%             
%             dS = dS + dS';
%             DS = -dS;
%             nele = length(eleobj);
%             for j=1:nele
%                 DS(j,j) = sum(dS(j,:))+dExt(j);
%             end
%             
%             exchS = DS*cS;
%             
%             exchM = D*cM;
%             
%             exchAct = D*conc;
%             
%             outdiff = (exchS-exchAct)*delta;
%             
%             cm = deltamud*cEM./(deltamud+ws.*sm);
%             dexp4  = sum(deltamud.*(conc-cEM))*delta;
%             bedMud = sum(M*(Vo-vf).*cb.*idx);
%             diff = dexp4 + bedMud;
%             
%             %upate Mass balance
%             MassBalance = MassBalance-diff;
        end
    end
end

