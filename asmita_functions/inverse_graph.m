function outgraph = inverse_graph(ingraph)
%
%-------function help------------------------------------------------------
% NAME
%   inverse_graph.m
% PURPOSE
%   To reverse the direction of a directed graph by taking the inverse of
%   the the adjacency matrix
% USAGE
%   outgraph = inverse_graph(ingraph)
% INPUTS
%   ingraph - digraph object that defines a flow field
% OUTPUTS
%   outgraph - digraph object that defines the reverse flow field
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    nn = numnodes(ingraph);
    [s,t] = findedge(ingraph);
    pthAdj = sparse(s,t,ingraph.Edges.Weight,nn,nn);
    outgraph = digraph(pthAdj');
    outgraph.Nodes.EleID = ingraph.Nodes.EleID;
    outgraph.Nodes.Type  = ingraph.Nodes.Type;
    outgraph.Nodes.Name  = ingraph.Nodes.Name;
    %reverse the in/out names 
    outgraph.Nodes.Name{end}  = '';  %add blank because must be unique
    outgraph.Nodes.Name{1}  = ingraph.Nodes.Name{end};
    outgraph.Nodes.Name{end}  = ingraph.Nodes.Name{1};
end