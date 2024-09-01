function outgraph = rescale_graph(ingraph,exchIn,isbalance)
%
%-------function help------------------------------------------------------
% NAME
%   rescale_graph.m
% PURPOSE
%   rescale the exchanges in a network based on the inputs defined in the
%   exchIn vector
% USAGE
%   outgraph = rescale_graph(ingraph,exchIn,isbalance)
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

    %get the matrix and exchanges for the input graph
    [exmatrix,oldIn,exOut,nodetxt] = graph2matrix(ingraph,nele);
    
    %update the exchange graph with the new inputs
    outgraph = matrix2graph(exmatrix,exchIn,exOut,nodetxt);
    
    if isbalance
        %impose a mass balance on each node of the new network (in=out). 
        %This assumes that the existing network has a mass balance at 
        %each node
        for j=1:2 %check outer and inner domains
            inid = find(exchIn(:,j)~=oldIn(:,j)); %inputs that have changed
            if all(inid) && sum(exchIn(:,j))==0, continue; end
            for i=1:length(inid)
                %for each input follow its flowpath updating the weights 
                %for all edges in the path
                idi = inid(i)+1; %id of the input node being updated
                                 %+1 accounts for the 'outside' node
                dFlow = exchIn(inid,2)-oldIn(inid,2);
                outgraph = update_successor_nodes(outgraph,dFlow,idi,idi);
            end
        end
    end
end
%%
function outgraph = update_successor_nodes(outgraph,dFlow,idnode,visited)
    %recursively check the nodes downstream and update the flows at each
    %node that is found
    outgraph = update_single_node(outgraph,dFlow,idnode);
    downstream = successors(outgraph,idnode); %nodes downstram of the current node
    if ~isempty(downstream)
        %call the functions recursively for the downstream nodes
        for i=1:length(downstream)
            idd = downstream(i);
            %ensure that the recursion does not get stuck in a loop
            visited = [visited,idd]; %#ok<AGROW> 
            if sum(visited==idd)<2
                %point has not already been updated
                outgraph = update_successor_nodes(outgraph,dFlow,idd,visited);                
            end            
        end
    end
end
%%
function outgraph = update_single_node(outgraph,dFlow,idnode)
    %update the weights based on a change in one of the input edges
    inp_ids = inedges(outgraph,idnode);
    out_ids = outedges(outgraph,idnode);
    inputs = outgraph.Edges.Weight(inp_ids); %new incoming flows
    outputs = outgraph.Edges.Weight(out_ids);%old outgoing flows
    intotal = sum(inputs);
    outotal = sum(outputs);
    if intotal~=outotal
        %flows do not balance so scale output by the change in total flow
        if length(outputs)>1
            %multiple outputs - scale change based on flow ratio of outputs
            for i=1:length(outputs)
                ratioFlow = outputs(i)/outotal;
                outgraph.Edges.Weight(out_ids(i)) = outputs(i)+dFlow*ratioFlow;
            end
        else
            %single output - pass through new value
            outgraph.Edges.Weight(out_ids) = intotal;
        end
    end
end
