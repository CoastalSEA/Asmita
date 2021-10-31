function flowgraph = matrix2graph(exmatrix,exchIn,exchOut,nodetxt,namemsg)
%
%-------function help------------------------------------------------------
% NAME
%   matrix2graph.m
% PURPOSE
%   use an exchange matrix and a vector of external exchanges to construct
%   a directed graph with labels defined in nodetxt
% USAGE
%   flowgraph = matrix2graph(exmatrix,exchange,nodetxt,namemsg)
% INPUTS
%   exmatrix - [nxn] matrix that defines the exchanges within the network
%   exchIn  - 2xn array of the exchanges with or from the outside, or [].
%   exchOut - 2xn array of the exchanges to the outside, or [].
%             NB: the network defines an exchange between one or two
%             domains. These external domains are typically spatially
%             distinct relative to the network being above/below,
%             north/south, left/right, seaward/landward, outer/inner, etc
%             The first column vector defines any exchange to/from the
%             outer domain and the second column vector defines
%             any exhanges to/from the inner domain; where, for example, 
%             outer might be the sea and inner the river(s).
%   nodetxt - struct with fields for:
%               nid - node id
%               ntype - node type
%               nname - node name
%   namemsg - node names must be unique. optional warning message used if 
%             duplicates names found
% OUTPUTS
%   flowgraph - handle to digraph object for the network
% NOTES
%   digraph can have exchanges in and/or out to define networks with inputs
%   only or outputs only. The unwanted exchange variable is passed as empty
% SEE ALSO
%   graph_matrix.m, used in Asmita Estuary and Advection classes
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    if nargin<5
        namemsg = [];                 %no namin error message supplied
    end
    
    nele = size(exmatrix,1);          %number of internal nodes in network
    userdata = zeros(nele+2,nele+2);  %full matrix with external exchanges
    
    %check the ingoing exchanges NB: for ingoing flows from the outer
    %domain only, it is necessary to use a 2 column input with zeros if the
    %second column
    if isempty(exchIn)                %no ingoing flows defined
        exchIn = zeros(nele,2);
    elseif size(exchIn,2)==1          %only one ingoing flow supplied
        %assume only the inner domain ingoing exchange has been provided        
        exchIn(:,2) = zeros(nele,1);  %avoid overwriting column 1
        exchIn = fliplr(exchIn);      %correct order so that outer domain 
    end                               %is in column 1
    
    %check the outgoing exchanges 
    if isempty(exchOut)               %no outgoing flows defined
        exchOut = zeros(nele,2);
    elseif size(exchOut,2)==1         %only one outgoing flow supplied
        %assume only the outer domain outgoing exchange has been provided
        exchOut(:,2) = zeros(nele,1);  
    end

    %exchanges in and out of the network
    userdata(1,2:end-1) = exchIn(:,1);    %outer domain ingoing
    userdata(2:end-1,1) = exchOut(:,1);   %outer domain outgoing
    userdata(end,2:end-1) = exchIn(:,2);  %inner domain ingoing
    userdata(2:end-1,end) = exchOut(:,2); %inner domain outgoing
    userdata(2:end-1,2:end-1) = exmatrix; %system network

    if sum(sum(userdata))==0
        userdata = 0;                     %remove matrix if no data
    end

    g = digraph(userdata);                %create graph from matrix
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
end     