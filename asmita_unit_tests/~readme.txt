The unit test files seek to reproduce results from previous versions. 

asmita_command_line_script.m illustrates how to run the model from the command prompt, change variables and plot some results

There are two unit testing  options
To test the current version of the model use:
	UserUnitTesting.m, AsmitaTest.m and runasmitamodel.m.
	This option runs the model with the model setup in the 'muiASM model files' folder and tests the results held in the 'muiASM TestData' folder
To test a model against the results from AsmitaOO v1.4 or v1.5:
	Run the required model so that there is at least on Case (uses last run case)
	Use UserUnitCFmatTesting.m and AsmitaCFmatfiles.m.
	This option compares the results from a muiASM model case from the model 'muiASM model files' folder with the results helld in the 'AsmitaOO Data files' folder.

Modifed 26/05/2023 to be installed and run using the Asmita App.

Comparison of mui and OO tests - 23/11/21
The parameter settings are not all the same in the testing version file and the version provided as a demo in the example folder.
#Humber 3EM 
fails against the OO version because there was an error in the application of the nodal cycle. The results agree if this is changed back.
#Amelander 
agrees with the results from the cfFortran code rather than the previous version of the Amelander test model
#Venice 
saltmarsh volumes are different to the previous version. Due to changes to tidal model 
#Inlet 
did not agree but there was a bug in the initial assignement of DQ because of an incorrect order of assignment of input variables related to conc. Once corrected the two versions agreed.
#Severn 
checked for volume data and concentration and all match. New model updates computation of tidal pumping so new output is different.

Test models for Humber 1EM, Southampton 4EM and Yangtze 30EM are new and test single element, interventions and hydrodynamics respectively.

