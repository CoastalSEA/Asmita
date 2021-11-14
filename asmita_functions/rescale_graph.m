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
%   exchIn - 2xn array of the exchanges with or from the outside.
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
    if nargin<3, isbalance = false; end
    
    nele = size(exchIn,1);
    if size(exchIn,2)==1              %same as matrix2graph assumption
        %assume only the inner domain ingoing exchange has been provided        
        exchIn(:,2) = zeros(nele,1);  %avoid overwriting column 1
        exchIn = fliplr(exchIn);      %correct order so that outer domain
    end                               %is in column 1

    %get the matrix and exchages for the input graph
    [exmatrix,oldIn,exOut,nodetxt] = graph2matrix(ingraph,nele);
    
    %update the exchange matrix with the new inputs
    outgraph = matrix2graph(exmatrix,exchIn,exOut,nodetxt);
    
    if isbalance
        %impose a mass balance on each node of the new network (in=out). 
        %This assumes that the existing network has a mass balance at 
        %each node
        for j=1:2 %check outer and inner domains
            inid = find(exchIn(:,j)~=oldIn(:,j)); %inputs that have changed
            if all(inid) && sum(exchIn(:,j))==0, continue; end
            for i=1:length(inid)
                %for each input follow its flowpath upating the weights for all
                %edges in the path
                idi = inid(i)+1; %id of the input node being updated
                                 %+1 accounts for the 'outside' node
                outgraph = update_successor_nodes(outgraph,idi);
            end
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
