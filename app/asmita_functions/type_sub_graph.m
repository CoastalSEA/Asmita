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
%   ntype - character vector, cell array of character vector, or 
%           numeric/logical indices of the type to use 
%           to sub-sample the network based on the Node.Type property
% OUTPUTS
%   outgraph - handle to digraph object for the 'ntype' sub-network
% SEE ALSO
%   used in Asmita Reach class
% EXAMPLE
%   outgraph = type_sub_graph(ingraph,{'','Channel); 
%   returns a sub graph wiht the end nodes and the 'Channel' type nodes. 
%   'end nodes' are typically inputs to the system eg from Sea or River.
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    if iscell(ntype) || isstring(ntype) || ischar(ntype)
        idx = ismatch(ingraph.Nodes.Type,ntype); %could be replaced by matches
    end
    outgraph = subgraph(ingraph,idx);
end