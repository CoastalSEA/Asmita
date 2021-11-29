%% ASMITA functions
% Summary of functions available in the _asmita_functions_ folder. Use the Matlab(TM)
% help function in the command window to get further details of each
% function.

%% Saltmarsh functions
% Functions that implement generic aspects of the saltmarsh model
% 

%%
% * *bioenhancedsettling.m*
% - Adjust the vertical exchange due to biological enhanced settling
% for the specified depth
%
% * *concovermarsh.m*
% - Calculate time averaged concentration over the marsh for a range of
% depths and return a dstable of depth, concentration and submergence time
% for an emulated spring-neap cycle.
%
% * *interpdepthload.m*
% - Iterate to find the equilibrium depth and sediment load, qm(s^-1)
%
% * *marsh_conc.m*
% - Calculate the concentration over a marsh flat using method proposed 
% by Krone(1987).
%
% * *morris_biocoeffs.m*
% - Calculate the biomass coefficients for Bc=aD+bD^2+c where D is
% the depth at peak biomass
%
% * *morris_eqdepth.m*
% - Solve the Morris equation for equilibrium depth

%% Network functions
% Functions used to manipulate networks and graphs of flow pathways

%%
% * *ele2node.m*
% - find graph node ids for specified element ids
% 
% * *graph2matrix.m*
% - extract the exchange matrix and vectors of external exchanges from
% a directed graph with labels defined in nodetxt
%
% * *inverse_graph.m*
% - To reverse the direction of a directed graph by taking the inverse of
% the the adjacency matrix
%
% * *matrix2graph.m*
% - use an exchange matrix and a vector of external exchanges to construct
% a directed graph with labels defined in nodetxt
%
% * *rescale_graph.m*
% - rescale the exchanges in a network based on the inputs defined in the
% exchIn vector
%
% * *setgraph.m*
% - update the Network or Flow tabs with directed graphs of network
%
% * *setmatrix.m*
% - generate UI to edit dispersion or advection matrix
%
% * *setnodetext.m*
% - generate UI to edit dispersion or advection matrix 
%
% * *type_sub_graph.m*
% - sub-sample a graph network based on node type
%

%% Other functions
% Other functions used in Asmita

%%
% * *readinputdata.m*
% -read advection timeseries data from file. Read format is file specific
% but can handle multiple columns with space or tab separators
%
% * *sealevelrise.m*
% - Function to compute sea level rise using linear or exponential rate of
% change
%
% * *userprismcoeffs.m*
% - parameters to define volume-prism relationships of the form:
% volume = alpha*prism^beta or volume = alpha*tidal range^beta

%% See Also
% The <matlab:doc('asmita_classes') additional classes> used in Asmita and
% the <matlab:open_manual manual>, which provides further details of setup and 
% configuration of the model.
