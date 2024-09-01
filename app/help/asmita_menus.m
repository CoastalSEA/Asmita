%% Menu Options
% Summary of the options available for each drop down menu.

%% File
% * *New*: clears any existing model (prompting to save if not already saved) and a popup dialog box prompts for Project name and Date (default is current date). 
% * *Open*: existing Asmita models are saved as *.mat files. User selects a model from dialog box.
% * *Save*: save a file that has already been saved.
% * *Save as*: save a file with a new or different name.
% * *Exit*: exit the program. The close window button has the same effect.

%% Tools
% * *Refresh*: updates Cases tab.
% * *Clear all > Project*: deletes the current project, including all Setup data and all Cases.
% * *Clear all > Figures*: deletes all results plot figures (useful if a large number of plots have been produced).
% * *Clear all > Cases*: deletes all Cases listed on the Cases tab but does not affect the model setup.

%% Project
% * *Project Info*: edit the Project name and Date
% * *Cases > Edit Description*: user selects a Case to edit the Case description.
% * *Cases > Edit Data Set*: initialises the Edit Data UI for editing data sets.
% * *Cases > Save*: user selects a data set to be saved from a list box of Cases and the is then prompted to name the file. The data are written to an Excel spreadsheet. 
% * *Cases > Delete*: user selects Case(s) to be deleted from a list box of Cases and results are then deleted (model setup is not changed).
% * *Cases > Reload*: user selects a Case to reload as the current parameter settings.
% * *Cases > View settings*: user selects a Case to display a table listing the parameters used for the selected Case. 
% * *Export/Import > Export*: user selects a Case class instance to export as a mat file.
% * *Export/Import > Import*: user selects an exported Case class instance (mat file) to be loaded.
%%
% *NB*: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the *Project>Cases>Edit Data Set* option 
% to make a selection and then use the ‘Copy to Clipboard’ button to paste 
% the selection to the clipboard.

%% Setup
% * *Estuary>System Parameters*: dialogue to set external concentration, convergence and wind parameters
% * *Estuary>Water Levels*: dialogue to define water levels, and forcing conditions such as rate of sea level rise and changes in tidal range.
% * *Estuary>Dispersion*: dialogue table to define the horizontal exchange
% matrix, which is used to determine the element connectivity.
% * *Estuary>Dynamic exchanges*: option to define time varying horizontal
% exchanges. Applied to all dispersion links defined as NaN in the
% horizontal exchange matrix (currently only an option for a single time
% series).
% * *Estuary>Model Constants*: A number of constants are used in the model. Generally, the default values are appropriate but these can be adjusted and saved with the project if required.
%
% * *Element>Define Elements*: when setting up a new model a listing of the different types of element allows the user to define the number of each type element to include in the model (elements can subsequently be added or deleted). Once a set of elements have been defined the user is prompted to edit the default names for each element, change the element type and choose to make the element non-erodible (for cases where the element can infill with sediment but if cannot erode.
% * *Element>Element Parameters*: definition of elements is on a continuous loop. Select Cancel from the element list box to finish defining element properties. 
% * *Element>Load Parameters*: rather than enter data one element at a time, this option allows the properties to be loaded from a file. 
% * *Element>Add*: select element type from list. The Define Element dialog box is then displayed allowing the user to edit the element name and whether it is erodible or not. Element properties and any Dispersion or Advection then need to be added using the individual dialog boxes for these components.
% * *Element>Delete*: select element to be deleted from the list of model elements.
%
% * *Saltmarsh>Species Parameters*: any number of saltmarsh species can be entered and a check is made to ensure that the number of values entered for minimum depth, maximum depth maximum biomass and species productivity is consistent with the number of species specified.
% * *Saltmarsh>Equilibrium Marsh Depth*: this utility allows the user to explore the response of the marsh depth to the defined saltmarsh parameters and to see how this response varies by adjusting the biomass production rates for each of the species included
% * *Saltmarsh>Biomass distribution*: displays a plot of how the biomass
% production of each species is distributed relative to tidal flat elevation
% * *Saltmarsh>Marsh-flat Animation*: shows the variation in marsh
% elevation over time starting from a bare tidal flat based on the current
% saltmarsh and water level settings
%
% * *Rivers>River Inputs*: this component can only be added once some channel elements have been setup in the model. If no rivers are defined, the dialog box requires the Channel id number, along with the flow rate and the concentration density of the rivers sediment load. 
% * *Rivers>River Advection*: dialogue table to define river flow advection connectivity.
% * *Rivers>River Timeseries*: If there is more than one river input, user is prompted to select a river input channel. The user is then prompted for a text file with time series data to define river flow rates (which optionally can also include concentrations). 
% * *Rivers>Delete Input*: select the source to be deleted from the list of channels with river inputs.
%
% * *Drift>Drift  Inputs*: this component can only be added once some channel elements have been setup in the model. If no rivers are defined, the dialog box requires the Channel id number, along with the flow rate and the concentration density of the rivers sediment load. 
% * *Drift>Drift Advection*: dialogue table to define drift advection connectivity.
% * *Drift>Drift Timeseries*: to add a drift times series first define a
% source input element. The user is then prompted for a text file with time series data to define drift rates as volumes of sediment per year (m3/year). This can either be a single input drift which is then passed through the network, based on the partitioning given by the initial advection matrix definition. Alternatively, a rate can be specified for selected links in the path. 
% * *Drift>Delete Input*: select the source to be deleted from the list of channels with river inputs.
%
% * *Interventions>Add or Edit*: select element from list of elements and then use table to add interventions. Selection of elements is in a continuous loop. Select Cancel from the element list box to finish adding interventions. 
% * *Interventions>Clear*: deletes interventions from all elements.
% * *Interventions>Load File*: interventions can also be loaded from a file. Select an element from the list of elements and then select the file to load from the file selection dialogue. 
% * *Interventions>Change Sign*: changes the sign of interventions defined
% for the selected element (volume and area).
%
% * *Hydraulics*: setup options when including idealised hydraulics in the model.
% This is only really needed for multi-element models of large systems with
% a large and variable river discharge. To run the option requires the
% CSTmodel App to be installed.
% * *Hydraulics>Additional Parameters*: dialogue to define the parameters
% needed to determine the along-channel hydraulics using the CST model.
% * *Hydraulics>Create Look-up*: creates a table of along-channel
% variations in water levels and discharge as a function of river discharge
% input.
% * *Hydraulics>Summary Plot*: provides a summary plot of along channel 
% water levels and velocities, allowing the range
% of variability as a function of river discharge to be investigated (this 
% can also be viewed on the _Tides>Hydraulics_ tab).
%
% * *Run Parameters>Time Step*: the time step parameters comprise the time step (years or parts of a year), the number of time steps and the intervals at which the model output should be saved. 
% * *Run Parameters>Conditions*: a number of components can be defined in the model setup and the user can then select whether, or not, to include them in a particular model run
% * *Run Parameters>Select Eq. Coefficients*: select the set of equilibrium coefficients to be used from the list provided. 
% * *Run Parameters>Edit Eq. Coefficients*: edit the selected coefficients
% for each Element Type.
% * *Run Parameters>Edit Eq. Coefficients List*: the coefficient sets are defined in userprismcoeffs.m. If a new set of coefficients is added, the name of the set can be added

%% Run
% * *Check Input*: checks mass balance of any advection flows included and that sufficient components have been defined for the model to run. This is not a full diagnostic check, so it does not ensure that the model output is valid – simply that the model will run.
% * *Run Model*: runs model, prompts for a Case description, which is added to the listing on the Cases tab.
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.

%% Help
% * *Help>Documentation*: access the online documentation for Asmita.
% * *Help>Manual*: access the manual pdf file.
% * *Help>Theory*: access papers that explain the background theory.

%% See Also
% The <matlab:asm_open_manual manual> provides further details of setup and 
% configuration of the model.