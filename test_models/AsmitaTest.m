classdef AsmitaTest <  matlab.unittest.TestCase
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
                        'Amelander_cfFortran.mat',...
                        'Venice 9EM.mat',...
                        'Severn 16EM.mat',...
                        'Pagham 3EM.mat',...
                        'Inlet 7EMts.mat'};       
    end
%%
    methods (Test)
        %asmitaTest is an instance of the class AsmitaTest
        %UserCase defines which test data set to load
        %testFile and testModel define the model to run for comparison
        %called from UserUnitTesting by selecting the relevant test case 
        
        function RunModelTest(asmitaTest,testModel)
            %run test using Matlab unittest function
            %import matlab.unittest.fixtures.PathFixture
            
            testdataset = load(asmitaTest.UserCase{testModel});            
            if isfield(testdataset,'sobj') && isa(testdataset.sobj,'Asmita')
                %new test data set file format (parts of sobj)
                lobj = testdataset.sobj;
                useCase = length(lobj.Cases.DataSets.AsmitaModel);
                dst = getDataset(lobj.Cases,useCase,1);
                %first 6 properties for all time steps, all elements 
                expSolution = dst.DataTable{:,1:6};
                isoldmodel = false;
            else
                %original test data set file format
                expSolution = testdataset.expSolution; %file with sub-set saved               
                isoldmodel = true;
            end   
            
            actSolution = runasmitamodel(asmitaTest.testFile{testModel},isoldmodel);
            
            if isoldmodel
                expSolution = reshape(expSolution,size(actSolution));
            end
            
            figure;
            jj = 1;
            plot(1:size(actSolution,1),actSolution(:,jj),...
                                1:size(actSolution,1),expSolution(:,jj));                            
            legend('Model Solution','Data Solution')
            
            asmitaTest.verifyEqual(actSolution,expSolution,'RelTol',0.01);
        end
    end
end

