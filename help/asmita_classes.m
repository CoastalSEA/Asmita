%% Asmita Classes
% The ASMITA model is built using the <matlab:doc('muitoolbox') muitoolbox>
% and a number of model specific classes.

%%
% The ASMITA app uses the following classes:
%%
% * *Advection*: handles all types of advective flow, including rivers, drift and tidal pumping
% * *ASM_model*: implements ASMinterface to allow functions in AMinterface to be overloaded
% * *ASM_Plots*: inherits muiPlots and customises plots for the ASMITA application
% * *ASM_PlotsUI*: implements the muiDataUI interface for a range of plot types
% * *ASMinterface*: basic functionality of the ASMITA model
% * *Asmita*: implements muiModelUI for the ASMITA application
% * *AsmitaModel*: handles the initialisation, running and posting of model results
% * *CSThydraulics*: handles the initialisation and calling of the CSTmodel App
% * *Drift*: setup and provide access to littoral drift sediment inputs
% * *Element*: setup and provide access to element parameters
% * *EqCoeffParams*: select and access the equilibrium coefficients to be used in the model (uses definitions in userprismcoeffs.m)
% * *Estuary*: parameters and methods for system, such as dispersion and linkages
% * *Interventions*: setup and provide access to imposed changes to system geometry (volume and area)
% * *Reach*: setup and provide access to parameters for reaches (a channel element plus linked elements)
% * *River*: setup and provide access to river flow/sediment inputs
% * *RunConditions*: define the specific conditions for a model run
% * *RunProperties*: define the model run time parameters
% * *Saltmarsh*: setup and access to saltmarsh parameters and functions to define behaviour 
% * *WaterLevels*: setup and access to water level definitions, tidal parameters and changes over time

%%
% See also the <matlab:doc('asmita_functions') additional functions> used in 
% Asmita and the <matlab:asm_open_manual manual>, which provides further details 
% of setup and configuration of the model.
