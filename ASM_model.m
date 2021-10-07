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
            obj = getClassObj(mobj,'Inputs','ASM_model');
            if isempty(obj)
                obj = ASM_model;       
            end
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
    end
   
end