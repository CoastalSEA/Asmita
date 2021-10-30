classdef ASM_model < ASMinterface                         
%
%-------class help---------------------------------------------------------
% NAME
%   ASM_model.m
% PURPOSE
%	Model computes inlet/estuary volume changes for a single time step. 
% 	Called recursively by AsmitaSim
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
        % eqScaling   %scaling of equilibirum relative to initial value
        % export      %summation of net export over run (for mass balance)
        % intervent   %summation of interventions over run (for mass balance)
        % dWLvolume   %summation of volume changes due to water level change
        % SedMbal     %sediment mass balance updated at each time step
        % WatMbal     %water mass balance updated at each time step
        % UniqueYears   %array of years for imposed changes
        % AnnualChange  %array of imposed volume and area changes for each element
        % DQ          %matix of dispersion and advection updated each time step
        % dqIn        %vector of 
    end
    
    methods (Access = private)
        function obj = ASM_model()                    
            %class constructor
        end
    end      
%%
    methods (Static, Sealed)                
        function setASM_model(mobj) %abstract method in ASMinterface
            %function to initialise an instance of ASM_model - called by
            %AsmitaSim.InitialiseModel
%             obj = getClassObj(mobj,'Inputs','ASM_model');
%             if isempty(obj)
%                 obj = ASM_model;       
%             end
            obj = ASM_model; 
            eleobj = getClassObj(mobj,'Inputs','Element');
            obj.eqScaling = ones(length(eleobj),1); %initialise default values
            setClassObj(mobj,'Inputs','ASM_model',obj);
        end
    end
%% ------------------------------------------------------------------------
    methods (Static)
    %overload functions in ASMinterface as required
 
        function asmitaEqFunctions(mobj)
            %function to define the equilibrium volume for a given prism
            eleobj = getClassObj(mobj,'Inputs','Element');
            ecpobj = getClassObj(mobj,'Inputs','EqCoeffProps');

            etypalpha = ecpobj.alpha;
            etypbeta = ecpobj.beta;
            eqType  = ecpobj.eqtype;
            %
            eletype = getEleProp(eleobj,'transEleType');
            %
            prism = Reach.getReachEleProp(mobj,'UpstreamPrism');
            HWL = Reach.getReachEleProp(mobj,'HWlevel');
            LWL = Reach.getReachEleProp(mobj,'LWlevel');
            %
            %equilibrium surface area taken as intial area (ie fixed)
            EqSA = getEleProp(eleobj,'InitialSurfaceArea');
            %Equilibrium depth over marsh elements
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
        function setDQmatrix(mobj,offset)
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
            [Q,qIn,~] = Advection.getAdvectionFlow(mobj,'River');
            [Qtp,qtpIn,~] = Advection.getAdvectionFlow(mobj,'Qtp');
            [Qs,qsIn,~] = Advection.getAdvectionFlow(mobj,'Drift');
            
            isQtpOnly = false;

            if isQtpOnly
                rncobj = getClassObj(mobj,'Inputs','RunConditions');
                if rncobj.IncTidalPumping
                    DQ = Qtp;
                    dEqIn = kCeI.*qtpIn;
                else
                    DQ = D+Q;
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
                
            %---------------------
            % may need to modify to use D+Q if river advection only (no tidal
            % pumping) and Qtp if tidal pumping included ie. QTP expresses
            % tidal and river exchange == D+Q
            %---------------------
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
            end

            setClassObj(mobj,'Inputs','ASM_model',obj);
        end        
        
    end
   
end