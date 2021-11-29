%% Tabs
% To examine what has been set-up the Tabs provide a summary of what is 
% currently defined. Note: these only update when clicked on using a mouse 
% and the values cannot be edited from the Tabs.

%% Tab options
% * *Cases*: lists the cases that have been run with a case id and description.
% * *Settings*: tabulates the current model settings (display only). 
% This tab has the following sub-tabs:
%%
% <html>
% <ul><ul>
% <li>System - tabulates the system and water level parameters.</li>
% <li>Elements - tabulates the key parameters of each element.</li>
% <li>Saltmarsh - tabulates the saltmarsh parameters and plots the species distribution biomass</li>
% <li>Interventions -  plots the volume and surface area changes defined for individual elements.</li>
% <li>Run Parameters - tabulates the key run parameters, conditions and equilibrium coefficients).</li>
% </ul></ul>
% </html>
%%
% * *Tides*: graphics of the tidal network, water levels and hydraulics (if included).
% This tab has the following sub-tabs:
%%
% <html>
% <ul><ul>
% <li>Network - graphic of the element network showing the connectivity and 
% horizontal exchanges between all elements (m3/s).</li>
% <li>Water Levels - plots the water levels over the run duration based on 
% the current settings.</li>
% <li>Hydraulics - plots the summary hydraulics (water levels and velocities)
% along channel if Hydraulic parameters have been defined.</li>
% <li>TP Discharge - plots the tidal pumping discharge as a function of distance from the mouth.</li>
% <li>TP Network - network graph plot of the tidal pumping advection
% flows.</li>
% </ul></ul>
% </html>
%%
% * *Flows*: graphics of the elements that have an advection (e.g. channels) 
% showing the connectivity and the flow rate (m3/s). This tab has the
% following sub-tabs:

%%
% <html>
% <ul><ul>
% <li>Rivers - network graph plot of the river advection flows.</li>
% <li>Drift - network graph plot of the drift advection flows.</li>
% <li>Input Summary - tabulated summary of the input nodes and the type of
% input.</li>
% <li>River Input - summary plot of any timeseries river advection inputs
% defined (flow and concentration if both are defined).</li>
% <li>Drift Input - summary plot of any timeseries drift advection inputs
% defined.</li>
% </ul></ul>
% </html>

%%
% * *Response*: summary of morphological response times based on system definition.
% * *Q-Plot*: user selects a variable from imported data or model output, and a plot is created on the tab. The plot cannot be edited.
% * *Stats*: tabulates the results from statistics UI.

%% Accessing Case meta-data
% On the *Cases* tab, using the mouse to click on a case record generates a
% table figure with a summary of the meta-data for the selected case. The
% figure includes buttons to Copy the data to the clipboard, view the
% DSproperties of the selected dataset and examine the Source of the
% dataset, which may be alist of files for imported data or details of the
% model used.

%% See Also
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model..