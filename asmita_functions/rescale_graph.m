function outgraph = rescale_graph(ingraph,exchIn,isbalance)
%
%-------function help------------------------------------------------------
% NAME
%   rescale_graph.m
% PURPOSE
%   rescale the exchanges in a network based on the inputs defined in the
%   exchIn vector
% USAGE
%   outgraph = rescale_graph(ingraph,exchIn)
% INPUTS
%   ingraph - an existing network graph that is to be rescaled
%   exchIn - new inputs to the network
%   isbalance - impose a mass balance at each node by scaling exchanges in
%               proportion to existing ratios
% OUTPUTS
%   outgraph - handle to digraph object for the rescaled network
% SEE ALSO
%   used in Asmita Advection class
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2021
%--------------------------------------------------------------------------
%
    outgraph = ingraph;
    %get the existing exchange matrix from the ingraph
    nn = numnodes(ingraph);
    [s,t] = findedge(ingraph);
    exchmatrix = sparse(s,t,ingraph.Edges.Weight,nn,nn);
    oldIn = full(exchmatrix(end,:));
    
    %update the exchange matrix with the new inputs
    if iscolumn(exchIn), exchIn = exchIn'; end   %force a row vector
    exchmatrix(end,2:end-1) = exchIn(2:end-1);
    g = digraph(exchmatrix);
    outgraph.Edges.Weight = g.Edges.Weight;
    
    if isbalance
        %impose a mass balance on each node of the new network (in=out). 
        %This assumes that the existing network has a mass balance at 
        %each node
        inid = find(exchIn~=oldIn); %inputs that have changed
        for i=1:length(inid)
            %for each input follow its flowpath upating the weights for all
            %edges in the path
            idi = inid(i); %id of the input node being updated
            outgraph = update_successor_nodes(outgraph,idi);
        end
    end
end
%%
function outgraph = update_successor_nodes(outgraph,idnode)
    %recursively check the nodes downstream and update the flows at each
    %node that is found
    outgraph = update_single_node(outgraph,idnode);
    downstream = successors(outgraph,idnode); %nodes downstram of the current node
    if ~isempty(downstream)
        %call the functions reqursively for the downstream nodes
        for i=1:length(downstream)
            idd = downstream(i);
            outgraph = update_successor_nodes(outgraph,idd);
        end
    end
end
%%
function outgraph = update_single_node(outgraph,idnode)
    %update the weights baed on a change in one of the input edges
    inp_ids = inedges(outgraph,idnode);
    out_ids = outedges(outgraph,idnode);
    inputs = outgraph.Edges.Weight(inp_ids); %new incoming flows
    outputs = outgraph.Edges.Weight(out_ids);%old outgoing flows
    intotal = sum(inputs);
    outotal = sum(outputs);
    if intotal~=outotal
        %flows do not balance so scale output by the change in total flow
        dFlow = intotal/outotal;
        for i=1:length(outputs)
            outgraph.Edges.Weight(out_ids(i)) = outputs(i)*dFlow;
        end
    end
end
