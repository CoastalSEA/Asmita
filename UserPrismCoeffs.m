function [alpha,beta,eqtype] = UserPrismCoeffs(UserSelection)
    % Parameters to define volume-prism relationships of the form:
    % volume = alpha*prism^beta.
    % eqtype flags whether to use prism (1) or tidal range (0)
    % UserSelection is selected from list using:
    % Setup>Run parameters>Equilibirum Coefficients
    % The list is defined in property UserEqCoeffOptions of RunProps
    %
    % Nov 2016: default list extended to include Flood/Ebb delta and Delta flat
    %
    %---------------------------------------------------------------------
    % AUTHOR
    % Ian Townend
    %
    % COPYRIGHT
    % CoastalSEA, (c) 2016
    %----------------------------------------------------------------------
    %  
    
    %define default structs and then overwrite with case specific values
    alpha = struct('Delta',1,'Channel',1,'Tidalflat',1,'Saltmarsh',1,...
                   'Storage',1,'FloodDelta',1,'EbbDelta',1,'DeltaFlat',1);
    beta = alpha;
    eqtype = alpha;

switch UserSelection       
    case 'Generic'                    % Generic parameter settings
        alpha.Delta = 5.0e-3;         %reworked W&A(1976) data
        alpha.Channel = 5.0e-3;       %COE paper Part 1 (Townend etal, 2016)
        alpha.Tidalflat = 0.08;       %COE paper Part 1 (Townend etal, 2016)
        alpha.Saltmarsh = 0.08;       %only used if there is no saltmarsh biomass
        alpha.Storage = 1;
        alpha.FloodDelta = 7000;      %Powell etal (2006)   
        alpha.EbbDelta = 5.0e-3;      %reworked W&A(1976) data
        alpha.DeltaFlat = 0.08;       %same as tidal flats
        %
        beta.Delta = 1.23;            %Walton & Adams (19760
        beta.Channel = 1.23;          %COE paper Part 1 (Townend etal, 2016)
        beta.Tidalflat = 1.08;        %COE paper Part 1 (Townend etal, 2016)
        beta.Saltmarsh = 1.08;        %only used if there is no saltmarsh biomass
        beta.Storage = 1;
        beta.FloodDelta = 0.37;       %Powell etal (2006) 
        beta.EbbDelta = 1.23;         %reworked W&A(1976) data
        beta.DeltaFlat = 1.08;        %same as tidal flats
    
    case 'Venice'                     % Venice 9EM model parameters
        alpha.Delta = 5.0e-3;
        alpha.Channel = 0.3;
        alpha.Tidalflat = 0.15;
        alpha.Saltmarsh = 0.3;
        %
        beta.Delta = 1.23;
        beta.Channel = 1;
        beta.Tidalflat = 1;
        beta.Saltmarsh = 1;
        
    case 'Amelander'                  % Wang's Amelander parameters
        alpha.Delta = 2.92157e-3;
        alpha.Channel = 10.241e-6;
        alpha.Tidalflat = 1.312e8;
        %
        beta.Delta = 1.23;
        beta.Channel = 1.55;
        beta.Tidalflat = 0;
        
     case 'Pagham'                    % Pagham parameters
        alpha.Delta = 8.5e-4;
        alpha.Channel = 1;
        alpha.Tidalflat = 1;
        alpha.Saltmarsh = 1;
        %
        beta.Delta = 1.23;
        beta.Channel = 1;
        beta.Tidalflat = 1;
        beta.Saltmarsh = 1;
end

