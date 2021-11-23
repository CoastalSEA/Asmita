The unit test files seek to reproduce resutls from previous versions. 
The parameter settings are not all the same in the testing version file and the version provided as a demo in the example folder.
Humber 3EM fails against the OO version because there was an error in the application of the nodal cycle. The results agree if this is changed back.
Amelander agrees with the results frim the cfFortran code rather than the previous version of the Amelander test model
Venice saltmarsh volumes are different to the previous version. Due to changes to tidal model 
Inlet did not agree but there was a bug in the initial assignement of DQ because of an incorrect order of assignment of input varialbes related to conc. Once corrected the two versions agreed.
Severn checked for volume data and concentration and all match. New model updates computation of tidal pumping so new output is different.
Test models for Humber 1EM, Southampton 4EM and Ynagtze 30EM are new and test single element, interventions and hydrodynamics respectively.