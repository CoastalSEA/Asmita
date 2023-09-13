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
        export      %summation of net export over run (for mass balance)
        import      %summation of net import over run (for mass balance)
        intervent   %summation of interventions over run (for mass balance)
        dWLvolume   %summation of volume changes due to water level change (for mass balance)
        SedMbal     %sediment mass balance updated at each time step
        WatMbal     %water mass balance updated at each time step
        UniqueYears   %array of years for imposed changes
        AnnualChange  %array of imposed volume and area changes for each element
        isFixed       %array of logical flags, true if change is non-erodible
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
            sa = getEleProp(eleobj,'SurfaceArea');
            n  = getEleProp(eleobj,'TransportCoeff');
            cb = getEleProp(eleobj,'BedConcentration');            
            ws = getEleProp(eleobj,'transVertExch');

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
            estobj  = getClassObj(mobj,'Inputs','Estuary');
            if rncobj.IncSaltmarsh
                Bk = Saltmarsh.BioProduction(mobj);
                dvb = diag(sign(n))*(Bk*vm)*robj.delta;
                if rncobj.IncDynamicElements && ws(estobj.ExchLinks(:,2))==0
                    %no biological production for elements with ws=0
                    dvb(estobj.ExchLinks(:,2)) = 0;
                end
                [dv_ero,ds_ero] = Saltmarsh.getEdgeErosion(mobj,robj);
            else
                dvb = zeros(size(dvf)); dv_ero = 0; ds_ero = 0;
            end
            %
            % if ~isempty(vis0)
            %     dvf(vis0) = 0; %force dvf to zero if vm=0
            % end
            
            % change in water volume - applies to both sediment and water
            % volumes hence use sign(n) rather than n>0 
            dvm = sign(n).*sa.*dwl;  
            if rncobj.IncDynamicElements && ws(estobj.ExchLinks(:,2))==0
                %no volume change for elements with ws=0
                dvm(estobj.ExchLinks(:,2)) = 0;
            end

            %
            % elenames = getEleProp(eleobj,'EleName');
            % table(dwl,conc,dvf,dvb,dvm,ve,vm,sm,Gam,dd,B,'RowNames',elenames)
            %
            vm = vm +dvm +dvf -dvb +dv_ero; %total change (moving surface)
            vf = vf +dvf -dvb +dv_ero;      %morphological change (fixed surface)
            vb = vb +dvb;                   %saltmarsh organic sedimentation
                                            %vm modifies bioproduction in
                                            %next time step            
            %check that elements have not infilled
            if any(vm(:)<=0)      %when elements go to zero retain small
                vf(vm<=0) = 999;  %value to prevent matrix becoming poorly 
                vm(vm<=0) = 999;  %conditioned     
            end       
            
            %update surace areas
            sa = ASM_model.asmitaAreaChange(mobj,robj,ds_ero);

            %depth values
            dm = vm./sa; df = vf./sa; de = ve./sa;
            
            %check that elements have not been removed (zero area)
            if any(sa(:)<=0)
                dm(sa<=0) = 0;
                df(sa<=0) = 0;
                de(sa<=0) = 0;
            end

            % assign results
            for i=1:nele
                eleobj(i).MovingVolume = vm(i);
                eleobj(i).FixedVolume = vf(i);
                eleobj(i).EleConcentration = conc(i);
                eleobj(i).BioProdVolume = vb(i);
                eleobj(i).MovingDepth = dm(i);
                eleobj(i).FixedDepth = df(i);
                eleobj(i).EqDepth = de(i);
            end
            %mass balance check
            ASM_model.MassBalance(mobj,robj);
            %diagnostic message to command window during runtime
            % t = robj.Time/mobj.Constants.y2s;
            % obj = getClassObj(mobj,'Inputs','ASM_model');
            % smb = obj.SedMbal; wmb = obj.WatMbal;
            % sprintf('t = %g; Sediment Balance %g; Water Balance %g',t,smb,wmb)
            setClassObj(mobj,'Inputs','Element',eleobj);
        end

%%
        function sa = asmitaAreaChange(mobj,robj,ds_ero)
            %calculate the changes in surface area over a timestep
            eleobj = getClassObj(mobj,'Inputs','Element');
            sa = getEleProp(eleobj,'SurfaceArea');
            se = getEleProp(eleobj,'EqSurfaceArea');

            %at the moment only a correction for marsh erosion is made

            erose = num2cell(sa+ds_ero);
            [eleobj.SurfaceArea] = erose{:};
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
%%
        function [B,dd] = BddMatrices(mobj)
            %compute the solution matrix B and vector dd (ie B*gamma-dd)
            %get model parameters
            eleobj = getClassObj(mobj,'Inputs','Element');
            sm = getEleProp(eleobj,'SurfaceArea');
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
            %****overloaded in ASM_model
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj  = getClassObj(mobj,'Inputs','Element');
            cE = getEleProp(eleobj,'EqConcentration');
            kCeI = River.getRiverProp(mobj,'tsRiverConc')./cE;
            %set up dispersion and advection matrices using dynamic
            %DispersionGraph and FlowGraph if IncDynamicElements is true
            %otherwise use Estuary-Dispersion and Advection-Flow properties
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
            %set vertical exchange model parameters
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            eleobj = getClassObj(mobj,'Inputs','Element');
            nele = length(eleobj);
            vm = getEleProp(eleobj,'MovingVolume');
            Vo = getEleProp(eleobj,'InitialVolume');
            ve = getEleProp(eleobj,'EqVolume');
            cb = getEleProp(eleobj,'BedConcentration');
            n  = getEleProp(eleobj,'TransportCoeff');

            %reset transVertExch to initial value because ws used to 
            %calculate B and dd to get unconstrained dvf
            [eleobj.transVertExch] = eleobj(:).VerticalExchange;

            %adjust ws if dynamic element with time varying exchange            
            estobj  = getClassObj(mobj,'Inputs','Estuary');
            if rncobj.IncDynamicElements
                tsyear = robj.DateTime/mobj.Constants.y2s;
                idx = find(estobj.DynamicExchange.Year<=tsyear,1,'last');             
                eleobj(estobj.ExchLinks(:,2)).transVertExch = ...
                                    estobj.DynamicExchange.Vertical(idx);  
            end  

            %adjust vertical exchange if saltmarsh included (function uses
            %transVertExch based on initial VerticalExchange settings 
            %modified if DynamicElements are active, ie the inorganic values)
            if rncobj.IncSaltmarsh
                wsbio = num2cell(Saltmarsh.BioSettlingRate(mobj));
                [eleobj.transVertExch] = wsbio{:};
            end

            %unadjusted morphological change in volume 
            [B,dd] = ASM_model.BddMatrices(mobj); %uses transVertExch
            Gam  = (ve./vm).^n;
            dvf  = (B*Gam-dd)./cb*robj.delta;
            
            %check if any erosion is greater than available sediment
            ws = getEleProp(eleobj,'transVertExch');
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
            ecpobj = getClassObj(mobj,'Inputs','EqCoeffParams');
            
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
            EqSA = getEleProp(eleobj,'SurfaceArea');
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
                isTReq = ~logical(eqType.(eletype{i})); %switch to true if tidal range equilbrium
                eleobj(i).EqSurfaceArea = EqSA(i);
                switch eletype{i}
                    case 'Saltmarsh'
                        if Deq(i)>0
                            eleobj(i).EqVolume = EqSA(i)*Deq(i);
                        else
                            eleobj(i).EqVolume = alpha*prism(i)^beta;
                        end
                    otherwise
                        if isTReq %appplies to any element type (eg tidalflat)
                            eleobj(i).EqVolume = alpha*(HWL(i)-LWL(i))^beta;
                        else
                            eleobj(i).EqVolume = alpha*prism(i)^beta;
                        end
                end
            end
            setClassObj(mobj,'Inputs','Element',eleobj);
        end
       
%%
        function conc = asmitaConcentrations(mobj,varargin)
            %get concentration taking acount of the horizontal exchanges
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj = getClassObj(mobj,'Inputs','Element');
            
            vm = getEleProp(eleobj,'MovingVolume');
            ve = getEleProp(eleobj,'EqVolume');
            sm = getEleProp(eleobj,'SurfaceArea');
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
             obj.import = 0;
             obj.intervent = 0;
             obj.dWLvolume = 0;
             obj.SedMbal = 0;
             obj.WatMbal = 0;
         end
         
%%
        function MassBalance(mobj,robj)
            obj = getClassObj(mobj,'Inputs','ASM_model');
            if robj.Time==0
                obj = ASMinterface.clearMassBalance(mobj);
            end
            % Load data
            delta = robj.delta;                    %time step in seconds
            time = robj.Time;                      %elapsed time in seconds

            eleobj = getClassObj(mobj,'Inputs','Element');
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
            sm = getEleProp(eleobj,'SurfaceArea');
            dwl = getEleProp(eleobj,'EleWLchange');
            [~,dExt] = Estuary.getDispersion(mobj); %uses ReachGraph
            [~,qIn,qOut] = Advection.getAdvectionFlow(mobj,'River'); %water flux in m3/s
            [~,qtpIn,qtpOut] = Advection.getAdvectionFlow(mobj,'Qtp'); %water flux in m3/s
            [~,qsIn,qsOut] = Advection.getAdvectionFlow(mobj,'Drift'); %sed.flow in m3/s

            % Calculate sediment mass balance (+ve is import to water column)
            impriv = sum(qIn.*cI)*(time>0);      %river import
            impqtp = sum(qtpIn.*cE)*(time>0);    %tidal pumping marine import
            impsed = sum(qsIn.*cE)*(time>0);     %littoral drift import
            dimp = (impriv+impqtp+impsed)*delta; %import to system
            obj.import = obj.import+dimp;        %export to marine
            dexp1  = sum(qOut.*conc)*(time>0);   %advective outflow
            dexp2  = sum(qtpOut.*conc)*(time>0); %tidal pumping outflow
            dexp3  = sum(qsOut.*conc)*(time>0);  %littoral drift export
            dexp4  = sum(dExt.*(conc-cE))*(time>0);%diffusion exchange with env (coarse)
            dexp   = (dexp1+dexp2+dexp3+dexp4)*delta;  %sum over a time step
            obj.export = obj.export+dexp;        %export to marine
            bed    = sum(M*(Vo-vf).*cb);         %bed changes
            obj.intervent = obj.intervent+sum(dV.*cb);  %user defined changes
            biosed = sum(vb.*cb);                %saltmarsh organic sedimentation
            %balance of cumulative changes from t=0 
            sedbal = obj.import+biosed-bed-obj.export-obj.intervent;
            obj.SedMbal = sedbal/sum(M*(Vo))*100;  %percentage change
            % obj.SedMbal = sedbal;

            % change in water volume - applies to both sediment and water
            obj.dWLvolume = obj.dWLvolume+sum((sm.*dwl));
            %change in water volume due to volume changes (water + sediment)
            dvm = sum(M*(vm-vf));
            estvol = dvm;
            
            %check balance
            obj.WatMbal = (estvol-obj.dWLvolume)/sum(M*(Vo))*100; %percentage change
            % obj.WatMbal = (estvol-obj.dWLvolume);

            setClassObj(mobj,'Inputs','ASM_model',obj);            
        end

    end
end

