function outgraph = type_sub_graph(ingraph,ntype)
%
%-------function help------------------------------------------------------
% NAME
%   type_sub_graph.m
% PURPOSE
%   sub-sample a graph network based on node type
% USAGE
%   outgraph = type_sub_graph(ingraph,ntype)
% INPUTS
%   ingraph - an existing network graph that is to be sub-sampled
%   ntype - character vector or numeric/ogical indices of the type to use 
%           to sub-sample the network based on the Node.Type property
% OUTPUTS
%   outgraph - handle to digraph object for the 'ntype' sub-network
% SEE ALSO
%   used in Asmita Reach class
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    if ischar(ntype)
        idx = strcmp(ingraph.Nodes.Type,ntype);
    end
    outgraph = subgraph(ingraph,idx);
end