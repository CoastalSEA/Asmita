classdef AsmitaTest <  matlab.unittest.TestCase
    % Test the running of Asmita against selected test cases
    % UserCase defines which test data set to load
    % testModel is used for the UI model selection
    % testFile defines the model to be run for comparison

    properties(TestParameter)        
        testModel = struct('H1EMTest',1,...
                    'HumberTest',2,...
                    'AmelanderTest',3,...
                    'VeniceTest',4,...
                    'SevernTest',5,...
                    'InletTest',6,...
                    'SouthamptonTest',7,...
                    'YangtzeTest',8);    
        UserCase = {'H1EMTestData',...
                    'HumberTestData',...
                    'AmelanderTestData',...
                    'VeniceTestData',...
                    'SevernTestData',...
                    'InletTestData',...
                    'SouthamptonTestData',...
                    'YangtzeTestData'};
        testFile = {'Humber 1EM.mat',...
                    'Humber 3EM.mat',...
                    'Amelander_cfFortran.mat',...
                    'Venice 9EM.mat',...
                    'Severn 16EM.mat',...
                    'Inlet 7EMts.mat',...
                    'Southampton 4EM.mat',...
                    'Yangtze 30EM allInts.mat'};  
    end
    %
    properties    
        ncol = 1;   %define which column of output array to plot
    end
%%
    methods (Test)
        %asmitaTest is an instance of the class AsmitaTest
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
                expSolution = dst.DataTable{:,:};
                isoldmodel = false;
            else
                %original test data set file format
                expSolution = testdataset.expSolution; %file with sub-set saved               
                isoldmodel = true;
            end   
            
            [actSolution,varlist,ele] = runasmitamodel(asmitaTest.testFile{testModel},...
                                                                isoldmodel);
            
            if isoldmodel
                %new model introduces dept variable so only use the first 3
                %variables for volumes (moving, fixed , equilibirum)
                expSolution = expSolution(:,:,1:3);
                expSolution = reshape(expSolution,size(actSolution));
            end
            
            %diagnostic plots             
            %Failure listing column subscript can be used to select jj
            jj = asmitaTest.ncol; 
            figure('Name','UnitTest','Tag','PlotFig');            
            plot(1:size(actSolution,1),actSolution(:,jj),...
                                1:size(actSolution,1),expSolution(:,jj));                            
            legend('Model Solution','Data Solution')
            var = varlist{ceil(jj/length(ele))};
            jjele = ele{rem(jj,length(ele))};
            title(sprintf('%s (%s)',var,jjele))
            
            asmitaTest.verifyEqual(actSolution,expSolution,'RelTol',0.01);
        end
    end
end

