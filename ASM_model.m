classdef ASM_model < ASMinterface                         
%
%-------class help---------------------------------------------------------
% NAME
%   ASM_model.m
% PURPOSE
%	Model computes inlet/estuary volume changes for a single time step. 
% 	Called by AsmitaModel
%  	The core model is provided by the methods in ASMinterface.
% NOTES
%	ASM_model inherits from ASMinterface allowing any of the core 
%  	functions to be overloaded. This enables research development of 
%  	alternative methods without interfering with core functionality.
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties (Transient) %inherited from ASMinterface:
        % export      %summation of net export over run (for mass balance)
        % import      %summation of net import over run (for mass balance)
        % intervent   %summation of interventions over run (for mass balance)
        % dWLvolume   %summation of volume changes due to water level change (for mass balance)
        % SedMbal     %sediment mass balance updated at each time step
        % WatMbal     %water mass balance updated at each time step
        % UniqueYears   %array of years for imposed changes
        % AnnualChange  %array of imposed volume and area changes for each element
        % isFixed       %array of logical flags, true if change is non-erodible
        % DQ          %matix of dispersion and advection updated each time step
        % dqIn        %vector of input dispersion and advection updated at each time step
        % dqOut       %vector of output dispersion and advection updated at each time step
    end
    
    methods (Access = private)
        function obj = ASM_model()                    
            %class constructor
        end
    end      
%%
    methods (Static, Sealed)                
        function setASM_model(mobj) %abstract method in ASMinterface
            %function to initialise an instance of ASM_model
            obj = ASM_model; 
            setClassObj(mobj,'Inputs','ASM_model',obj);
        end
    end
%% ------------------------------------------------------------------------
    methods (Static)
    %overload functions in ASMinterface as required 
        function asmitaEqFunctions(mobj)
            %function to define the equilibrium volume for a given prism,
            %tidal range or drift and equilbirum surface area if variable
            eleobj = getClassObj(mobj,'Inputs','Element');
            rncobj  = getClassObj(mobj,'Inputs','RunConditions');
            inp = ASM_model.eqInputParams(mobj,eleobj);

            %Equilibrium depth over marsh elements
            if rncobj.IncSaltmarsh
                inp.Deq = Saltmarsh.EqDepth(mobj); 
            else
                %marsh element present but biology not active so use
                %specified equilibrium volume for marsh flat.
                inp.Deq = -ones(size(inp.eletype));
            end  
            
            %set equilibrium volumes and areas fixed surface area uses 
            %prism or tidal range power laws
            ASM_model.asmitaEqVolume(mobj,eleobj,rncobj,inp);
        end

%%          
        function asmitaEqVolume(mobj,eleobj,rncobj,inp)
            %function to define the equilibrium volume for a given prism,
            %tidal range or drift. Equilbirum surface area is initialised 
            %in Element.initialiseElements and used here as a constant 
            %value (inp.EqSA).

            %assign equilibrium volume to each element 
            for i=1:length(eleobj)
                alpha = inp.alpha.(inp.eletype{i});
                beta = inp.beta.(inp.eletype{i});
                tidalrange = inp.TR(i);
                isTReq = ~logical(inp.eqType.(inp.eletype{i})); %switch to true if tidal range equilibrium
                switch inp.eletype{i}
                    case 'Saltmarsh'
                        if inp.Deq(i)>0
                            %depth within species range
                            EqVol = inp.EqSA(i)*inp.Deq(i);
                        elseif inp.Deq(i)==0
                            %no water depth over marsh
                            EqVol = 0;
                        elseif inp.Deq(i)==-1
                            %depth greater than maximum species depth, 
                            %or root to Morris equation not found. NB: in 
                            %this formulation the eq.coefficients need to 
                            %be adjusted to suit site. 
                            if isTReq %appplies to any element type (eg tidalflat)
                                EqVol = alpha*inp.EqSA(i)*tidalrange^beta;
                            else
                                EqVol = alpha*inp.UPrism(i)^beta;
                            end                            
                        end

                    case {'Beachface','Shoreface','Spit','DeltaFlat'}
                        %scale beach equilbrium as function of drift rate, 
                        %tidal range or prism
                        EqVol = ASM_model.beachEqVolume(mobj,eleobj,rncobj,inp,i);
                    otherwise
                        if isTReq %appplies to any element type (eg tidalflat)
                            EqVol = alpha*inp.EqSA(i)*tidalrange^beta;
                        else
                            EqVol = alpha*inp.UPrism(i)^beta;
                        end
                end
                if EqVol<0, EqVol=0; end
                eleobj(i).EqVolume = EqVol;
            end
            setClassObj(mobj,'Inputs','Element',eleobj);
        end 

%%
        function EqVol = beachEqVolume(mobj,eleobj,rncobj,inp,ii)
            %set the beach, shoreface, delta flat equilibrium volumes
            %set beach equilbrium as function of drift rate, 
            %or tidal range or prism
            alpha = inp.alpha.(inp.eletype{ii});
            beta = inp.beta.(inp.eletype{ii});
            tidalrange = inp.TR(ii);
            isTReq = ~logical(inp.eqType.(inp.eletype{ii})); %switch to true if tidal range equilibrium
            if rncobj.IncDrift && rncobj.IncDriftTS && ~isTReq
                %for time varying drift rate use changes in drift rate
                advobj = getClassObj(mobj,'Inputs','Advection');
                [DR,initGraph] = getFlowRatio(advobj,mobj,'Drift'); %DR is struct of ratio (q/q0), diff (difference, q-q0) and diffratio (difference/initial value - dq/q0)                           
                elename = getEleProp(eleobj,'EleName');
                DiffRatio = DR.diffratio(strcmp(initGraph.Nodes.Name,elename(ii)));
                eleVol = getEleProp(eleobj,'InitialVolume');
                EqVol = eleVol(ii)*(1+alpha.*DiffRatio^beta);
            elseif isTReq
                %for constant drift rate use tidal range if eqtype=0
                EqVol = alpha*inp.EqSA(ii)*tidalrange^beta;
            else
                %for constant drift rate use tidal prism if eqtype=1
                EqVol = alpha*inp.UPrism(ii)^beta;
            end                        
        end

%%
        function inp = eqInputParams(mobj,eleobj)
            %get parameters used as input for the equilbirum volume and area
            ecpobj = getClassObj(mobj,'Inputs','EqCoeffParams');
        
            inp.alpha = ecpobj.alpha;
            inp.beta = ecpobj.beta;
            inp.eqType  = ecpobj.eqtype;
            %
            inp.eletype = getEleProp(eleobj,'transEleType');
            inp.EqSA = getEleProp(eleobj,'EqSurfaceArea');
            inp.UPrism = Reach.getReachEleProp(mobj,'UpstreamPrism');
            inp.RPrism = Reach.getReachEleProp(mobj,'ReachPrism');
            inp.TR = Reach.getReachEleProp(mobj,'TidalRange');
        end

%% 
        function ok = setDQmatrix(mobj,offset)
            % overloads version in ASMinterface
            % used to test D+Q+Qtp and just Qtp when tidal pimping included
            
            %set the DQ, dqIn and conc properties to be used for a time step
            %set up dispersion and advection matrices using dynamic
            %ReachGraph and FlowGraph if IncDynamicElements is true
            %otherwise use Estuary-Dispersion and Advection-Flow properties
            obj = getClassObj(mobj,'Inputs','ASM_model');
            eleobj  = getClassObj(mobj,'Inputs','Element');
            cE = getEleProp(eleobj,'EqConcentration');
            kCeI = River.getRiverProp(mobj,'tsRiverConc')./cE;
            [D,dExt] = Estuary.getDispersion(mobj);
            [Q,qIn,qOut] = Advection.getAdvectionFlow(mobj,'River');
            [Qtp,qtpIn,qtpOut] = Advection.getAdvectionFlow(mobj,'Qtp');
            [Qs,qsIn,qsOut] = Advection.getAdvectionFlow(mobj,'Drift');
            if isnan(Qs), ok = 0; return; end

            isQtpOnly = false;

            if isQtpOnly
                 %--------------------- NOT CURRENTLY WORKING -------------
                 %TEST code to examine use of Qtp in place of exchange
                 %components D+Q that relate to the channel elements
                 % *does not work for Severn but does work for Venice*
                 %---------------------------------------------------------
                rncobj = getClassObj(mobj,'Inputs','RunConditions');
                if rncobj.IncTidalPumping
                    DQ = D; dEqIn = dExt;
                    %needs checking - surely Qtp can be +/-? so need an
                    %index for non-zero values??
                    idx = abs(Qtp)>0;        %index of qtp exchanges
                    DQ(idx) = Qtp(idx);
                    idi = abs(qtpIn)>0;      %index of qtpIn exchanges
                    dEqIn(idi) = qtpIn(idi); %assumes no input from rivers 
                else                         %which would be negative qtpOut
                    DQ = D+Q;                %and need correcting by kCeI
                    dEqIn = dExt+kCeI.*qIn;
                end
                %
                switch offset
                    case 'flow+drift'    %include flow and drift
                        obj.DQ = DQ+Qs;
                        obj.dqIn = dEqIn+qsIn;
                    case 'flow'          %include flow only
                        obj.DQ = DQ;
                        obj.dqIn = dEqIn;
                    case 'drift'         %include drift only
                        obj.DQ = D+Qs;
                        obj.dqIn = dExt+qsIn;
                    otherwise            %use no advections (offset = 'none')
                        obj.DQ = D;
                        obj.dqIn = dExt;
                end                 
                
           
            else
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
            end

            setClassObj(mobj,'Inputs','ASM_model',obj);
            ok = 1;
        end        
    end
   
end