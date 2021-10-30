function test_graph(funcname,option)
%
%-------function help------------------------------------------------------
% NAME
%   test_graph.m
% PURPOSE
%   functions to test utility functions
% USAGE
%   test_graph(funcname)
% INPUT
%   funcname - name of function as a character vector or string
%   option - load an example of a advection or a dispersion matrix
% OUTPUT
%   runs test code for requested function
% 
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%    
    %load a test data set - mat file containing an instance of the ASMITA 
    %Advection  class
%     [sfile,spath] = uigetfile({'*.mat','MAT-files (*.mat)'},'Open testdata mat file');
%     if sfile==0, return; end
%     testdata = load([spath,sfile],'adv');
    spath = 'D:\Work\Tools\MATLAB\MUImodels2\muiApps\Asmita\asmita_functions\';
    switch option
        case 1
            sfile = 'severn_advection.mat';
            testdata = load([spath,sfile],'adv');
        case 2
            sfile = 'severn_dispersion.mat';
            testdata = load([spath,sfile],'disp');
    end
    
    %create input digraph
    ingraph = get_graph(testdata);
    
    switch funcname
        case 'exchange_graph'
            test_exchange(ingraph);
        case 'graph_matrix'
            test_matrix(ingraph,testdata);
        case 'inverse_graph'
            test_inverse(ingraph);        
        case 'rescale_graph'
            test_rescale(ingraph,exchIn,exchOut);
        case 'type_sub_graph'
            test_type(ingraph);
    end
end

%%
function test_exchange(ingraph)
    %test the exchange_graph function
    figure;
    plot(ingraph,'EdgeLabel',ingraph.Edges.Weight);   
end

%%
function test_matrix(ingraph,testdata)
    %test the extraction of defining matrix and vectors from a graph    
%     [exmatrix,exchIn,exchOut,nodetxt] = graph_matrix(ingraph); %#ok<ASGLU>
    
    if isfield(testdata,'adv')
        
        tdata1 = testdata.adv.matrix;
        tdata2 = testdata.adv.in;
        nele = size(tdata1,1);
        [exmatrix,exchIn,exchOut,nodetxt] = graph_matrix(ingraph,nele);
        %this is a submatrix
        
        exIn = exchIn(:,2);
	    exOut = exchOut(:,1);
    else
        [exmatrix,exchIn,exchOut,nodetxt] = graph_matrix(ingraph); 
        tdata1 = testdata.disp.matrix;
        tdata2 = testdata.disp.ext;
        exIn = exchIn(:,1);
        exOut = exchOut(:,1);
    end
    % tt = table(exmatrix,exchIn,exchOut)
    % check1 = sum(sum(tdata1-exmatrix))
    % check2 = sum(tdata2-exIn)
    outgraph = exchange_graph(exmatrix,exIn,exOut,nodetxt);
    plot_graph(ingraph,outgraph,'graph matrix');
end

%%
function test_inverse(ingraph)
    %test the inverse_graph function - reverse direction of a directed graph
    outgraph = inverse_graph(ingraph);
    plot_graph(ingraph,outgraph,'inverse graph');
end

%%
function test_rescale(ingraph,exchIn,exchOut)
    %test the rescale_graph function 
    idx = exchIn>0 | exchOut>0;
    exIn = [0;exchIn(idx);0];
    exIn = exIn*1.2;  %increase all inputs by 20%
    % exIn(4) = exIn(4)*1.2;   %increase only one input
    %test the application of new inputs without checking mass balance
    outgraph = rescale_graph(ingraph,exIn,false);    
    plot_graph(ingraph,outgraph,'rescale graph: no mass balance');
    %test altering inputs and imposing a mass balance
    outgraph = rescale_graph(ingraph,exIn,true);    
    plot_graph(ingraph,outgraph,'rescale graph: including balance');
end

%%
function test_type(ingraph)
    %test the subsampling of the network based on type
    ingraph.Nodes.Type(2:7) = {'Channel'};
    outgraph = type_sub_graph(ingraph,'Channel');
    plot_graph(ingraph,outgraph,'sub graph for Channel type');
end

%%
function ingraph = get_graph(testdata)
    %create a digraph based on the testdata
    if isfield(testdata,'adv')
        exmatrix = testdata.adv.matrix;
        exchIn = testdata.adv.in;
        exchOut = testdata.adv.out;
        names = {'Out','In'};
    else
        exmatrix = testdata.disp.matrix;
        exchIn = testdata.disp.ext;
        exchOut = [];
        names = {'','Ext'};
    end
    nodetxt = struct('nid',[],'ntype','','nname','');
    nele = size(exmatrix,1);
    %node inputs are padded at either end for the input and output
    nodetxt.nid = [0,1:nele,0]';   %must be a column vector
    nodetxt.ntype = repmat({'test'},nele+2,1);   %pad for ends
    nodetxt.nname = [names(1);repmat({'channel'},nele,1);names(2)]; %add ends
    ingraph = exchange_graph(exmatrix,exchIn,exchOut,nodetxt);
end

%%
function plot_graph(g1,g2,funcname)
    %plot before and after graphs
    figure;
    subplot(1,2,1)
    plot(g1,'EdgeLabel',g1.Edges.Weight);
    subplot(1,2,2)
    plot(g2,'EdgeLabel',g2.Edges.Weight);
    sgtitle(funcname)
end