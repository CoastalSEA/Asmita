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
%   option - load an example of (1) a advection or (2) a dispersion matrix
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
        case 'matrix2graph'
            test_exchange(ingraph);
        case 'graph2matrix'
            test_matrix(ingraph,testdata);
        case 'inverse_graph'
            test_inverse(ingraph);        
        case 'rescale_graph'
            test_rescale(ingraph,testdata);
        case 'type_sub_graph'
            test_type(ingraph);
        case 'ele2node'
            test_ele2node(ingraph);
    end
end

%%
function test_exchange(ingraph)
    %test the exchange_graph function
    figure;
    plot(ingraph,'EdgeLabel',ingraph.Edges.Weight);
    title('matrix2graph')
end

%%
function test_matrix(ingraph,testdata)
    %test the extraction of defining matrix and vectors from a graph    
    if isfield(testdata,'adv')
        tdata = testdata.adv.matrix;
        nele = size(tdata,1);        
        [exmatrix,exchIn,exchOut,nodetxt] = graph2matrix(ingraph,nele);
    else
        [exmatrix,exchIn,exchOut,nodetxt] = graph2matrix(ingraph); 
    end
    % tt = table(exmatrix,exchIn,exchOut)
    outgraph = matrix2graph(exmatrix,exchIn,exchOut,nodetxt);
    plot_graph(ingraph,outgraph,'graph2matrix');
end

%%
function test_inverse(ingraph)
    %test the inverse_graph function - reverse direction of a directed graph
    outgraph = inverse_graph(ingraph);
    plot_graph(ingraph,outgraph,'inverse graph');
end

%%
function test_rescale(ingraph,testdata)
    %test the rescale_graph function 
    if isfield(testdata,'adv')
        exchIn = testdata.adv.in;
    else
        exchIn = testdata.disp.ext;
    end
    exIn = exchIn*1.2;  %increase all inputs by 20%
    % exIn = exchIn;   exIn(4) = exchIn(4)*1.2;   %increase only one input
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
function test_ele2node(ingraph)
    %test extraction of node IDs from element IDs or Type
    figure; plot(ingraph);
    nodeids = ele2node(ingraph,[0,4,5,6]);
    figure; plot(subgraph(ingraph,nodeids));
    nodeids = ele2node(ingraph,{'Channel','test'});
    figure; plot(subgraph(ingraph,nodeids));
end
%%
function ingraph = get_graph(testdata)
    %create a digraph based on the testdata
    if isfield(testdata,'adv')
        exmatrix = testdata.adv.matrix;
        exchIn = testdata.adv.in;
        exchOut = testdata.adv.out;
        names = {'Out','In'};
        % exmatrix(12,3) = 44;  %used to test non-adjacent inputs in matrix
        % exchIn(12) = 44;      %however does not update mass-balance
    else
        exmatrix = testdata.disp.matrix;
        exchIn = testdata.disp.ext;
        exchIn(:,2) = zeros(length(exchIn),1);
        exchOut = [];
        names = {'Ext',''};
    end
    nodetxt = struct('nid',[],'ntype','','nname','');
    nele = size(exmatrix,1);
    %node inputs are padded at either end for the input and output
    nodetxt.nid = [0,1:nele,0]';   %must be a column vector
    nodetxt.ntype = [{'none'};repmat({'test'},nele,1);{'none'}];   %pad for ends
    nodetxt.nname = [names(1);repmat({'channel'},nele,1);names(2)];%add ends
    ingraph = matrix2graph(exmatrix,exchIn,exchOut,nodetxt);
end

%%
function plot_graph(g1,g2,funcname)
    %plot before and after graphs
    figure('Name','Graph plot','Tag','PlotFig');
    subplot(1,2,1)
    plot(g1,'EdgeLabel',g1.Edges.Weight);
    subplot(1,2,2)
    plot(g2,'EdgeLabel',g2.Edges.Weight);
    sgtitle(funcname)
end