function nodeid = ele2node(agraph,eleid)
%
%-------function help------------------------------------------------------
% NAME
%   ele2node.m
% PURPOSE
%   find graph node ids for specified element ids
% USAGE
%   nodeid = ele2node(agraph,eleid)
% INPUTS
%   agraph - the graph or digraph to use to select nodes
%   eleid - element IDs (numerical), or element type (character vector,
%           or cell array)
% OUTPUTS
%   nodeid = nodes in the graph with the corresponding eleid
% SEE ALSO
%   used in Asmita Reach class
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    if ischar(eleid) || iscell(eleid) || isstring(eleid)
        nodetype = agraph.Nodes.Type;
        nodeid = find(ismatch(nodetype,eleid));
    else
        nodelist = agraph.Nodes.EleID;
        %find the eleids that are members in nodelist
        [~,nid] = ismember(eleid,nodelist);
        nodeid = nid(nid>0); %remove any null members (0)
    end
end