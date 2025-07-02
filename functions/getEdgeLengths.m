% Returns the lengths of the edges of graph described by nodes and edges

% INPUTS
% nodes: a d x n array of vertex coordinates where d is the dimension and n
% is the number of vertices
% edges: a m x 2 array listing the indices of the vertices connected by
% each edge

% OUTPUTS
% lengths: a 1 x m array of the lengths of the edges

function lengths = getEdgeLengths(nodes, edges)
    nEdges = size(edges,1);
    lengths = nan(1,nEdges);
    
    for i = 1:nEdges        
        p1 = nodes(:,edges(i,1));
        p2 = nodes(:,edges(i,2));
        lengths(i) = norm(p1 - p2);
    end
end