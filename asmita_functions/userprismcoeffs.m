function [alpha,beta,eqtype] = userprismcoeffs(UserSelection)
%
%-------function help------------------------------------------------------
% NAME
%   userprismcoeffs.m
% PURPOSE
%   Parameters to define volume-prism relationships of the form:
%   volume = alpha*prism^beta or volume = alpha*tidal range^beta
% USAGE
%   [alpha,beta,eqtype] = userprismcoeffs(UserSelection)
% INPUTS
%   UserDelection - 
% OUTPUTS
%   The 3 outputs return structs containing values for each element type
%   alpha - scale parameter for prism relationship
%   beta - shape paraemter for prism relationship
%   eqtype - flags whether to use prism (1) or tidal range (0)
% NOTES
%   Setting the scale parameter, beta, to 0, has the effect of applying a
%   constant value defined by alpha.
%   The Element Types are as defined by the GeoTypes property in the Asmita
%   class and MUST be the same for selection to work. 

% SEE ALSO
%   EqCoeffParams class in Asmita handles the selection of a parameter set
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%

%define default structs and then overwrite with case specific values
alpha = struct('Channel',1,'Tidalflat',1,'Saltmarsh',1,...
               'Storage',1,'FloodDelta',1,'EbbDelta',1,'DeltaFlat',1,...
               'Beachface',0,'Shoreface',0,'Spit',0);
beta = alpha;
eqtype = alpha;

switch UserSelection       
    case 'Default'
        alpha.EbbDelta = 2.92157e-3;
        beta.EbbDelta = 1.23;
        
    case 'Generic'                    % Generic parameter settings
        alpha.Channel = 5.0e-3;       %COE paper Part 1 (Townend etal, 2016)
        alpha.Tidalflat = 0.08;       %COE paper Part 1 (Townend etal, 2016)
        alpha.Saltmarsh = 0.08;       %only used if there is no saltmarsh biomass
        alpha.Storage = 1;
        alpha.FloodDelta = 7000;      %Powell etal (2006)   
        alpha.EbbDelta = 5.0e-3;      %reworked W&A(1976) data
        alpha.DeltaFlat = 1;    
        alpha.Beachface = 1;
        alpha.Shoreface = 1;
        alpha.Spit = 1;
        %
        beta.Channel = 1.23;          %COE paper Part 1 (Townend etal, 2016)
        beta.Tidalflat = 1.08;        %COE paper Part 1 (Townend etal, 2016)
        beta.Saltmarsh = 1.08;        %only used if there is no saltmarsh biomass
        beta.Storage = 1;
        beta.FloodDelta = 0.37;       %Powell etal (2006) 
        beta.EbbDelta = 1.23;         %reworked W&A(1976) data
        beta.DeltaFlat = 1;   
        beta.Beachface = 1;
        beta.Shoreface = 1;
        beta.Spit = 1;
    
    case 'Venice'                     % Venice 9EM model parameters
        alpha.EbbDelta = 5.0e-3;
        alpha.Channel = 0.3;
        alpha.Tidalflat = 0.15;
        alpha.Saltmarsh = 0.3;
        %
        beta.EbbDelta = 1.23;
        beta.Channel = 1;
        beta.Tidalflat = 1;
        beta.Saltmarsh = 1;
        
    case 'Amelander'                  % Wang's Amelander parameters
        alpha.EbbDelta = 2.92157e-3;
        alpha.Channel = 10.241e-6;
        alpha.Tidalflat = 1.312e8;
        %
        beta.EbbDelta = 1.23;
        beta.Channel = 1.55;
        beta.Tidalflat = 0;
        
     case 'Pagham'                    % Pagham parameters
        alpha.EbbDelta = 8.5e-4;
        alpha.Channel = 1;
        alpha.Tidalflat = 1;
        alpha.Saltmarsh = 1;
        %
        beta.EbbDelta = 1.23;
        beta.Channel = 1;
        beta.Tidalflat = 1;
        beta.Saltmarsh = 1;
end

