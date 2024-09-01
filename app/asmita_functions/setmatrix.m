function [exmatrix,exchIn,exchOut,flowgraph] = setmatrix(eleobj,figtitle,...
                                                promptxt,inoutxt,inpmatrix)
%
%-------function help------------------------------------------------------
% NAME
%   setmatrix.m
% PURPOSE
%   generate UI to edit dispersion or advection matrix 
% USAGE
%   [exmatrix,exchIn,exchOut] = setmatrix(eleobj,figtitle,promptxt,inoutxt,inpmatrix)
% INPUT
%   eleobj - instance of Element class - used to define ID and name
%   figtitle - title of the UI figure
%   promptxt - descriptive text that preceeds table used to prompt user
%   inoutxt - cell array of labels for input and output (source and sink)
%   inpmatrix - data to be used to populate the table
% OUTPUT
%   exmatrix - [nxn] matrix that defines the exchanges within the network
%   exchIn  - 2xn array of the exchanges with or from the outside.
%   exchOut - 2xn array of the exchanges to the outside.
%             NB: the network defines an exchange between one or two
%             domains. These external domains are typically spatially
%             distinct relative to the network being above/below,
%             north/south, left/right, seaward/landward, outer/inner, etc
%             The first column vector defines any exchange to/from the
%             outer domain and the second column vector defines
%             any exhanges to/from the inner domain; where, for example, 
%             outer might be the sea and inner the river(s).
%   flowgraph - handle to digraph object for the network
% SEE ALSO
%   matrixtableUI.m and tablefigureUI.m
%
% Author: Ian Townend
% CoastalSEA (c)Oct 2021
%--------------------------------------------------------------------------
%
    nele = size(inpmatrix,1)-2;
    nodetxt = setnodetext(eleobj,inoutxt);
    %UI for user to edit exchanges
    newtable = matrixtableUI(figtitle,promptxt,nodetxt.nname,inpmatrix);
    g = digraph(newtable{:,:});
    %remove elements that are not connected
    [in,out] = findedge(g);
    idx = unique([in,out]);
    flowgraph = subgraph(g,idx);          %flowgraph for given exchanges
    %add id and names to nodes
    flowgraph.Nodes.EleID = nodetxt.nid(idx);  
    flowgraph.Nodes.Type = nodetxt.ntype(idx);
    %node names must be unique
    nname = check_unique_names(nodetxt.nname(idx),true,namemsg);
    flowgraph.Nodes.Name = nname;    
    %resultant exchange matrix and input/output arrays
    [exmatrix,exchIn,exchOut] = graph2matrix(flowgraph,nele);
end
%%
function msg = namemsg()
    msg1 = 'Duplicate element names are not allowed';
    msg2 = 'Some names have been modified';
    msg3 = 'Edit Element names using Setup>Elements>Define Elements';
    msg = sprintf('%s\n%s\n%s',msg1,msg2,msg3);
end