function nodetxt = setnodetext(eleobj,inoutxt)
%
%-------function help------------------------------------------------------
% NAME
%   setnodetext.m
% PURPOSE
%   generate UI to edit dispersion or advection matrix 
% USAGE
%   nodetxt = setnodetext(eleobj,inoutxt)
% INPUT
%   eleobj - instance of Element class (which need not be current settings)
%   inoutxt - cell array of labels for input and output (source and sink)
% OUTPUT
%   nodetxt - struct with fields for:
%               nid - node id
%               ntype - node type
%               nname - node name
% SEE ALSO
%   setmatrix.m, graph2matrix.m, matrix2graph.m
%
% Author: Ian Townend
% CoastalSEA (c)Oct 2021
%--------------------------------------------------------------------------
%
    if isempty(eleobj(1).transEleType)
        eletype = getEleProp(eleobj,'EleType');  
    else
        eletype = getEleProp(eleobj,'transEleType');  
    end
    eleid = getEleProp(eleobj,'EleID');
    elename = getEleProp(eleobj,'EleName');

    %use the element id, type and name to return the nodetxt, a struct with
    %nid - node id, ntype - node type, nname - node name
    nele = length(eleid);
    nodetxt.nid = zeros(nele+2,1);
    nodetxt.nid = [0;eleid;0];
    nodetxt.ntype = [{''};eletype;{''}];                   
    nodetxt.nname = [inoutxt(1);elename;inoutxt(2)];
end