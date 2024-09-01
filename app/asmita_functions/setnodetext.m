function nodetxt = setnodetext(eleobj,inoutxt,atype)
%
%-------function help------------------------------------------------------
% NAME
%   setnodetext.m
% PURPOSE
%   generate UI to edit dispersion or advection matrix 
% USAGE
%   nodetxt = setnodetext(eleobj,inoutxt,atype)
% INPUT
%   eleobj - instance of Element class (which need not be current settings)
%   inoutxt - cell array of labels for input and output (source and sink)
%   atype - cell array of types to be included (optional)
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
    
    if nargin>2 && ~isempty(atype)
        idx = ismatch(eletype,atype); %could be replaced by matches
        nele = sum(idx);
    else
        nele = length(eleid);
        idx = 1:nele;
    end

    %use the element id, type and name to return the nodetxt, a struct with
    %nid - node id, ntype - node type, nname - node name    
    nodetxt.nid = zeros(nele+2,1);
    nodetxt.nid = [0;eleid(idx);0];
    nodetxt.ntype = [{''};eletype(idx);{''}];                   
    nodetxt.nname = [inoutxt(1);elename(idx);inoutxt(2)];
end