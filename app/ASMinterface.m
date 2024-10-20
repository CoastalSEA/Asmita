classdef ASMinterface < handle
    % Abstract class for ASM_model. Has all functions needed to run model
    % all of which can be overloaded in ASM_model
    %
    %----------------------------------------------------------------------
    % NOTES
    % 1) logical array of 'illcondition' tests for matrix stability. If any
    %    test fails the run is terminated. The array comprises:
    %    1=DQ; 2=Sm; 3=W; 4=DQ+Sm*W
    %    ie horizontal exchanges, surface areas, vertical exchanges, and 
    %    coefficient matrix used in inverse operations. These tests are run
    %    for volume, surface area, and concentration calculations. Area tests
    %    1 and 3 use D/h, Q/h and W/h, hence depths may also be cause of illcondition.
    % 2) Surface areas that are <=0 are set to a value of 999 to
    %    avoid matrices becoming illconditioned.
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
        dqOut       %vector of output dispersion and advection updated at each time step
        errCode = 0 %see notes above.
    end

    methods (Abstract,Static)
        setASM_model(mobj)
    end
    
    methods (Static)
        %[alpha,beta,eqtype] = userprismcoeffs(UserSelection) %external function

%%
        function ok = asmitaVolumeChange(mobj,robj)
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
            se = getEleProp(eleobj,'EqSurfaceArea');
            n  = getEleProp(eleobj,'TransportCoeff');
            cb = getEleProp(eleobj,'BedConcentration');            

            %get concentration before updating volumes
            [conc,ok] = ASM_model.asmitaConcentrations(mobj);
            if ok<1, return; end

            [B,dd,ok] = ASM_model.BddMatrices(mobj);
            if ok<1, return; end
            Gam  = (ve./vm).^n;

            % morphological change in volume            
            dvf  = (B*Gam-dd)./cb*robj.delta;

            %biomass change in volume
            rncobj = getClassObj(mobj,'Inputs','RunConditions');
            if rncobj.IncSaltmarsh
                Bk = Saltmarsh.BioProduction(mobj);
                dvb = diag(sign(n))*(Bk*vm)*robj.delta;
                [dv_ero,~] = Saltmarsh.setEdgeErosion(mobj,robj);
            else
                dvb = zeros(size(dvf)); dv_ero = 0; %ds_ero = 0; 
            end
            
            % change in water volume - applies to both sediment and water
            % volumes hence use sign(n) rather than n>0 
            dvm = sign(n).*sa.*dwl;  

            % elenames = getEleProp(eleobj,'EleName');
            % table(dwl,conc,dvf,dvb,dvm,ve,vm,sm,Gam,dd,B,'RowNames',elenames)
            
            vm = vm +dvm +dvf -dvb +dv_ero; %total change (moving surface)
            vf = vf +dvf -dvb +dv_ero;      %morphological change (fixed surface)
            vb = vb +dvb;                   %saltmarsh organic sedimentation
                                            %vm modifies bioproduction in
                                            %next time step            
            %check that elements have not infilled
            if any(vm<=0)
                idx = find(vm<=0);
                msgtxt = sprintf('Volume of element %d <= 0\n',idx);
                warndlg(sprintf('%s\nAbort run in ASM_model.asmitaVolumeChange',msgtxt))
                ok = 0; return;
                %vf(vm<999) = 999;  %value to prevent matrix becoming poorly 
                %vm(vm<999) = 999;  %conditioned     
            end       

            %update surace areas (marsh edge erosion and variable area)
            % [sa,ok] = ASM_model.asmitaAreaChange(mobj,robj);  
            % if ok<1, return; end

            %depth values
            dm = vm./sa; df = vf./sa; de = ve./se;
            
            %check that elements have not been removed (zero area)
            %or eroded (sa<=999 set to 0)
            if any(sa<=999)
                idx = sa<=999;
                dm(idx) = 0;
                df(idx) = 0;
                de(idx) = 0;
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
            setClassObj(mobj,'Inputs','Element',eleobj);

            %mass balance check
            ok = ASM_model.MassBalance(mobj,robj);
            %diagnostic message to command window during runtime
            % t = robj.Time/mobj.Constants.y2s;
            % obj = getClassObj(mobj,'Inputs','ASM_model');
            % smb = obj.SedMbal; wmb = obj.WatMbal;
            % sprintf('t = %g; Sediment Balance %g; Water Balance %g',t,smb,wmb)            
        end

%%
        function [B,dd,ok] = BddMatrices(mobj)
            %compute the solution matrix B and vector dd (ie B*gamma-dd)
            %get model parameters
            obj  = getClassObj(mobj,'Inputs','ASM_model');
            eleobj = getClassObj(mobj,'Inputs','Element');
            sm = getEleProp(eleobj,'SurfaceArea');
            ws = getEleProp(eleobj,'transVertExch');
            n = getEleProp(eleobj,'TransportCoeff');
            cE = getEleProp(eleobj,'EqConcentration');
            
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
            illCondition = ASMinterface.checkMatrix(DQ,Sm,W);
            if any(isnan(DQ),'all') || any(isinf(DQ),'all')
                errordlg(sprintf('DQ matrix contains NaN or Inf in BddMatrices\nRun aborted'))
                B = 0;   dd = 0;  ok = 0;  
            elseif ~any(illCondition) || (illCondition(3) && ~illCondition(4)) 
                %if W is illconditioned but DQ+SW is not do not abort. This can
                %be caused by setting vertical exchange to zero when
                %simulating setback schemes           
                B    = CE*M*Sm*W*(I-(DQ+Sm*W)\Sm*W);
                dd   = CE*M*Sm*W*((DQ+Sm*W)\dqIn);
                ok = 1;
            else
                if obj.errCode==0, obj.errCode = illCondition; end                
                errtxt = num2cell(illCondition);
                errtxt = sprintf('DQ: %d; S: %d; W: %d; DQ+SW: %d',errtxt{:});
                msgtxt = sprintf('Volume matrix is illconditioned.\nError code: %s\nRun aborted',errtxt);
                errordlg(msgtxt)
                B = 0;   dd = 0;  ok = 0;                               
            end
            %full equation is as follows(DQ and dqIn replace terms in brackets)
            % B   = cE*M*Sm*Wz*(I-(D+Q+Qs+Qtp+Sm*Wz)\Sm*Wz);
            % dd  = cE*M*Sm*Wz*((D+Q+Qs+Qtp+Sm*Wz)\(dExt+qtpIn+kCeI.*qIn+qsIn));          
        end
        
%% 
        function ok = setDQmatrix(mobj,offset)
            %set the DQ, dqIn and conc properties to be used for a time step
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj  = getClassObj(mobj,'Inputs','Element');
            cE = getEleProp(eleobj,'EqConcentration');
            kCeI = River.getRiverProp(mobj,'tsRiverConc')./cE;
            %set up dispersion and advection matrices using dynamic
            %DispersionGraph and FlowGraph if IncDynamicElements is true
            %otherwise use Estuary-Dispersion and Advection-Flow properties
            [D,dExt] = Estuary.getDispersion(mobj);
            [Q,qIn,qOut] = Advection.getAdvectionFlow(mobj,'River');
            [Qtp,qtpIn,qtpOut] = Advection.getAdvectionFlow(mobj,'Qtp');
            [Qs,qsIn,qsOut] = Advection.getAdvectionFlow(mobj,'Drift');
            if isnan(Qs), ok = 0; return; end

            %to set eqCorV in Element.setEleAdvOffsets need to only
            %include some of the advections based on conditions set
            %and offset = RunConditions.Adv2Offset
            %whereas runtime calls use all advections included in run
            %and offset = RunConditions.Adv2Inc
            switch offset
                    case 'flow+drift'    %include flow and drift
                        obj.DQ = D+Q+Qtp+Qs;
                        obj.dqIn = dExt+kCeI.*qIn+qtpIn+qsIn;
                        obj.dqOut = dExt+qOut+qtpOut+qsOut;
                    case 'flow'          %include flow only
                        obj.DQ = D+Q+Qtp;
                        obj.dqIn = dExt+kCeI.*qIn+qtpIn;
                        obj.dqOut = dExt+qOut+qtpOut;
                    case 'drift'         %include drift only
                        obj.DQ = D+Qs;
                        obj.dqIn = dExt+qsIn;
                        obj.dqOut = dExt+qOut;
                    otherwise            %use no advections (offset = 'none')
                        obj.DQ = D;
                        obj.dqIn = dExt;
                        obj.dqOut = dExt;
            end   

            setClassObj(mobj,'Inputs','ASM_model',obj);
            ok = 1;
        end

%%
        function ok = setVertExch(mobj,robj)
            %set vertical exchange model parameters
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            eleobj = getClassObj(mobj,'Inputs','Element');
            sa = getEleProp(eleobj,'SurfaceArea');
            %reset transVertExch to initial value because ws used to 
            %calculate B and dd to get unconstrained dvf
            [eleobj.transVertExch] = eleobj(:).VerticalExchange;

            %adjust ws if dynamic element with time varying exchange            
            estobj  = getClassObj(mobj,'Inputs','Estuary');
            if rncobj.IncDynamicElements
                tsyear = robj.DateTime/mobj.Constants.y2s;
                idx = find(estobj.DynamicExchange.Year<=tsyear,1,'last'); 
                nrow = size(estobj.ExchLinks,1);
                %can be more than one NaN link but all use the same
                %exchange table (ie same changes area made to all NaN links)
                news = num2cell(repmat(estobj.DynamicExchange.Vertical(idx),1,nrow));
                [eleobj(estobj.ExchLinks(:,2)).transVertExch] = news{:}; 
                %update instance for use in BioSettlingRate
                setClassObj(mobj,'Inputs','Element',eleobj);
            end  

            %adjust vertical exchange if saltmarsh included (function uses
            %transVertExch based on initial VerticalExchange settings 
            %modified if DynamicElements are active, ie the inorganic values) 
            if rncobj.IncSaltmarsh
                wsbio = Saltmarsh.BioSettlingRate(mobj);
                wsbio = num2cell(wsbio);
                [eleobj.transVertExch] = wsbio{:};
                %update instance for use in BddMatrices
                setClassObj(mobj,'Inputs','Element',eleobj);
            end

            %check for inactive elements
            if any(sa==999)
                ids = sa==999;
                saws = num2cell(ones(1,sum(ids))*eps(100));
                [eleobj(ids).transVertExch] = saws{:};
            end
            %check if any erosion is greater than available sediment for
            %any non-erodible elements
            ok = ASM_model.checkFixity(mobj,robj);
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
            TR = Reach.getReachEleProp(mobj,'TidalRange');
            %
            %equilibrium surface area taken as intial area (ie fixed)
            EqSArea = getEleProp(eleobj,'SurfaceArea');
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
                EqSA = EqSArea(i);
                switch eletype{i}
                    case 'Saltmarsh'
                        if Deq(i)>0
                            EqVol = EqSA*Deq(i);
                        else
                            EqVol = alpha*prism(i)^beta;
                        end
                    otherwise
                        if isTReq %appplies to any element type (eg tidalflat)
                            EqVol = alpha*TR(i)^beta;
                        else
                            EqVol = alpha*prism(i)^beta;
                        end
                end
            end
            if EqVol<0, EqVol=0; end
            eleobj(i).EqVolume = EqVol;
            setClassObj(mobj,'Inputs','Element',eleobj);
        end

%%
        function [conc,ok] = asmitaConcentrations(mobj,varargin)
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
            illCondition = ASMinterface.checkMatrix(DQ,Sm,W);
            if any(isinf(Gam)) || any(isnan(Gam))
                %error in Gam - divide by zero if any(vm=0)
                errordlg(sprintf('Some volumes are zero in asmitaConcentrations\nRun aborted'))
                conc = zeros(size(cE));                  ok = 0;  
            elseif ~any(illCondition) || (illCondition(3) && ~illCondition(4)) 
                %if W is illconditions but SW is not do not abort. This can
                %be caused by setting vertical exchange to zero when
                %simulating setback schemes   
                conc = (DQ+Sm*W)\(Sm*W*Gam+(dqIn)).*cE;  ok = 1;                
            else
                if obj.errCode==0, obj.errCode = illCondition; end
                errtxt = num2cell(illCondition);
                errtxt = sprintf('DQ: %d; S: %d; W: %d; DQ+SW: %d',errtxt{:});
                msgtxt = sprintf('Concentration matrix is illconditioned.\nError code: %s\nRun aborted',errtxt);
                errordlg(msgtxt)
                conc = zeros(size(cE));                  ok = 0;  
            end
            %once DQ and dqIn are expanded, this is equivalent to:
            % cnm = (D+Q+Qs+Qtp+Sm*Wz)\(Sm*Wz*Gam+(dE+qtpin+kCeI.*qEin+qSin))*cE;
        end

%%
        function ok = checkFixity(mobj,robj)
            % check whether any element is fixed and if so check that 
            %change does not exceed sediment availabel relative to initial 
            %volume taking account of any fixed interventions        
            eleobj = getClassObj(mobj,'Inputs','Element');
            nele = length(eleobj);
            vm = getEleProp(eleobj,'MovingVolume');
            Vo = getEleProp(eleobj,'InitialVolume');
            ve = getEleProp(eleobj,'EqVolume');
            cb = getEleProp(eleobj,'BedConcentration');
            n  = getEleProp(eleobj,'TransportCoeff');

            %unadjusted morphological change in volume 
            [B,dd,ok] = ASM_model.BddMatrices(mobj); %uses transVertExch
            if ok<1, return; end
            Gam  = (ve./vm).^n;
            dvf  = (B*Gam-dd)./cb*robj.delta;            
            %check if any erosion is greater than available sediment for
            %any element
            ws = getEleProp(eleobj,'transVertExch');
            vf = getEleProp(eleobj,'FixedVolume');   
            vfixint = getEleProp(eleobj,'eqFixedInts');
            noero = ~getEleProp(eleobj,'Erodible'); %ie if Erodible false 
            %current sediment availability relative to reference volume
            %SedAvailable sign(n).*(Vref-vf)<0 - element eroding
            Vref = Vo+vfixint(:,1);
            SedAvailable = sign(n).*(Vref-vf); 
            %when noero true and dvf>0 ie erodiing
            %correction=1 if SedAvailable>dvf
            %correction=0 if SedAvailable<=0
            %correction=SedAvailable/SedChange if SedAvailable<=dvf
            EleChange = sign(n).*dvf; %+ve -> more water or less sediment
            for j=1:nele
                correction = 1;                           %SedAvailable>dvf 
                if noero(j) && EleChange(j)>0 
                    if SedAvailable(j)<=0                 %SedAvailable<=0
                        correction = eps(100);   
                    elseif SedAvailable(j)<EleChange(j)   %SedAvailable<dvf
                        correction = SedAvailable(j)/EleChange(j);    
                    end
                end
                %update vertical exchange based on any correction needed
                eleobj(j).transVertExch = ws(j)*correction;                
            end
            setClassObj(mobj,'Inputs','Element',eleobj);
            ok = 1;            
        end

%%
        function illcondition = checkMatrix(amat1,amat2,amat3)
            %check whether matrices are ill conditioned for inversion. Return logical 
            %array for check of each matrix and combination amat1+amat2*amat3. 
            %use 'fro' because Frobenius norm is better for sparse matrices
            threshold = 1/eps;
            amatrix = amat1+amat2*amat3;
        
            illcondition = [cond(amat1,'fro')>threshold,...
                            cond(amat2,'fro')>threshold,...
                            cond(amat3,'fro')>threshold,...
                            cond(amatrix,'fro')>threshold];
        end

%%      
        function obj = clearMassBalance(mobj)
            %clear properties used to compute mass balance
            obj = getClassObj(mobj,'Inputs','ASM_model');
            obj.export = 0;
            obj.import = 0;
            obj.intervent = 0;
            obj.dWLvolume = 0;
            obj.SedMbal = 0;
            obj.WatMbal = 0;
        end
         
%%
        function ok = MassBalance(mobj,robj)
            obj = getClassObj(mobj,'Inputs','ASM_model');
            if robj.Time==0
                obj = ASMinterface.clearMassBalance(mobj);
            end
            % Load data
            delta = robj.delta;                    %time step in seconds
            time = robj.Time;                      %elapsed time in seconds

            eleobj = getClassObj(mobj,'Inputs','Element');
            cE = getEleProp(eleobj,'EqConcentration');
            % cES = estobj.EqConcCoarse;
            % cEM = estobj.EqConcFine;           
            cb = getEleProp(eleobj,'BedConcentration');
            cI = River.getRiverProp(mobj,'tsRiverConc');
            conc = getEleProp(eleobj,'EleConcentration');
            n = getEleProp(eleobj,'TransportCoeff');
            M = diag(sign(n));
            wM = diag(n>0); %#ok<NASGU> 
            
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
            [Qs,qsIn,qsOut] = Advection.getAdvectionFlow(mobj,'Drift'); %sed.flow in m3/s
            if isnan(Qs), ok = 0; return; end

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
            obj.intervent = obj.intervent+sum(dV.*cb);
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
            ok = 1;
        end

    end
end

