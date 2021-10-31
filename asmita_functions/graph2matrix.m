function [exmatrix,exchIn,exchOut,nodetxt] = graph2matrix(flowgraph,nele)
%
%-------function help------------------------------------------------------
% NAME
%   graph2matrix.m
% PURPOSE
%   extract the exchange matrix and vectors of external exchanges from
%   a directed graph with labels defined in nodetxt
% USAGE
%   [exmatrix,exchIn,exchOut,nodetxt] = graph2matrix(flowgraph)
% INPUTS
%   flowgraph - handle to digraph object for the network
%   nele - number of nodes required in the matrix being returned (optional)
%          if not specified the number of nodes in the flowgraph is used 
%          to determine nele.
% OUTPUTS
%   exmatrix - [nxn] matrix that defines the exchanges within the network
%   exchIn  - 2xn vector of the exchanges with or from the outside.
%   exchOut - 2xn vector of the exchanges to the outside.
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
% NOTES
%   digraph can have exchanges in and/or out to define networks with inputs
%   only or outputs only. The unwanted exchange variable is passed as empty
% SEE ALSO
%   exchange_graph.m and rescale_graph.m, used in Asmita Estuary and 
%   Advection classes
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%

    %recover the node text data
    nodetxt.nid = flowgraph.Nodes.EleID;
    nodetxt.ntype = flowgraph.Nodes.Type;
    nodetxt.nname = flowgraph.Nodes.Name;
    
    if nargin>1 && ~isempty(nele)    %user has specified the network size
        nn = nele;
    else
        nn = numnodes(flowgraph);    %use the flowgraph to determine size
    end
    
    [s,t] = findedge(flowgraph);
    spmatrix = sparse(s,t,flowgraph.Edges.Weight,nn,nn);
    fullmatrix = full(spmatrix);     %matrix of network and exchanges
    
    idx = nodetxt.nid>0;             %~idx can can be 1, end, or [1,end]
    exmatrix = fullmatrix(idx,idx);  %exchanges within the network
    
    if sum(~idx)==1
        %when only input or output included, the matrix is not symmetric
        %pad the exchanges to return [nx2] array
        exchIn = zeros(size(exmatrix,1),2); exchOut = exchIn;
        exchIn(:,~idx) = fullmatrix(~idx,idx)';  %exchanges into the network
        exchOut(:,~idx) = fullmatrix(idx,~idx);  %exchanges out of the network
        %add the missing exchange to the nodetxt to maintain symmetry
        ids = find(~idx);
        if ids==1                    %outer exists so add inner            
            nodetxt.nid(end+1) = 0;
            nodetxt.ntype{end+1} = '';
            nodetxt.nname{end+1} = '';
        else                         %inner exists so add outer
            nodetxt.nid = [0;nodetxt.nid];
            nodetxt.ntype = [{''};nodetxt.ntype];
            nodetxt.nname = [{''};nodetxt.nname];            
        end
    else
        exchIn = fullmatrix(~idx,idx)';  %exchanges into the network
        exchOut = fullmatrix(idx,~idx);  %exchanges out of the network
    end

    
    matrixsze = size(exmatrix,1);
    %if the matrix is a subset of the network expand to full matrix
    %based on the number of elements required equal to nele
    if nargin>1 && ~isempty(nele) && nele>matrixsze
        %initialise arrays
        exmatrix = zeros(nele,nele);                     %null matrix
        exchIn = zeros(nele,2); exchOut = exchIn;        %null vectors
        newid = zeros(nele,1); newtype = repmat({''},nele,1); newname = newtype;
        %assign values based on node ID stored in nodetxt.nid
        nids = nodetxt.nid(idx);
        exmatrix(nids,nids) = fullmatrix(idx,idx);
        exchIn(nids,:) = fullmatrix(~idx,idx)';
        exchOut(nids,:) = fullmatrix(idx,~idx);
        %update nodetext to the new array dimension
        newid(nids) = nodetxt.nid(idx);
        newtype(nids) = nodetxt.ntype(idx);
        newname(nids) = nodetxt.nname(idx); 
        nodetxt.nid = [nodetxt.nid(1);newid;nodetxt.nid(end)];
        nodetxt.ntype = [nodetxt.ntype(1);newtype;nodetxt.ntype(end)];
        nodetxt.nname = [nodetxt.nname(1);newname;nodetxt.nname(end)];
    end
end