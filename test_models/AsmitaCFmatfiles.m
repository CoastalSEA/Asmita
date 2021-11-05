classdef AsmitaCFmatfiles < matlab.unittest.TestCase
    % Test the running of Asmita against selected test cases
    properties(TestParameter)
            UserCase = {'HumberTestData',...
                        'AmelanderTestData',...
                        'VeniceTestData',...
                        'SevernTestData',...
                        'PaghamTestData',...
                        'InletTestData'};
            testModel = struct('HumberTest',1,...
                        'AmelanderTest',2,...
                        'VeniceTest',3,...
                        'SevernTest',4,...
                        'PaghamTest',5,...
                        'InletTest',6);
            testFile = {'Humber 3EM.mat',...
                        'Amelander 3EM.mat',...
                        'Venice 9EM.mat',...
                        'Severn 16EM.mat',...
                        'Pagham 3EM.mat',...
                        'Inlet 7EM.mat'};       
    end
%%
    methods (Test)
        %asmitaCFmatfiles is an instance of the class AsmitCFmatfiles
        %UserCase defines which test data set to load
        %testFile defines the results to be used in the comparison
        %the original AsmitaOO files have the same name and are held in a
        %subfoler: AsmitaOO Data files.
        %called from UserUnitTesting by selecting the relevant test case 
        
        function RunModelTest(asmitaCFmatfiles,testModel)
            %compare with actual soluation from a file saved for the model
            testdataset = load(asmitaCFmatfiles.testFile{testModel});
            lobj = testdataset.sobj;
            useCase = length(lobj.Cases.DataSets.AsmitaModel);
            dst = getDataset(lobj.Cases,useCase,1);
            %first 6 properties for all time steps, all elements 
            actSolution = dst.DataTable{:,1:6};
            %need to remove totals column from each variable
            ncol = size(actSolution,2);
            ntot = ncol/6;
            nint = ntot:ntot:ncol;
            actSolution(:,nint) = [];
            
            clear useCase dst lobj testdataset
            
            %load exepcted results from TestData file
            testdataset = load([pwd,'\AsmitaOO Data files\',asmitaCFmatfiles.UserCase{testModel},'.mat']);
            %first 3 properties for all time steps, all elements 
            expSolution = testdataset.expSolution; %file with sub-set saved
            expSolution = reshape(expSolution,size(actSolution));
            clear testdataset

            figure;
            plot(1:size(actSolution,1),expSolution(:,1));
            legend('Model Solution','Data Solution')
            asmitaCFmatfiles.verifyEqual(actSolution,expSolution,'RelTol',0.05);
        end
    end
end