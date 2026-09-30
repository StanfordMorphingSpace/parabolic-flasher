function adj_faces = getAdjacentFaces(vert_u, edges, faces)

adj_face = cell(length(edges), 1);
% rearrange how faces are stored
for i = 1:length(edges)
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);

        % Bending forces
        faces_with_node1        = (faces(:, 1) == p1_index) | (faces(:, 2) == p1_index) | (faces(:, 3) == p1_index);
        faces_with_node2        = (faces(:, 1) == p2_index) | (faces(:, 2) == p2_index) | (faces(:, 3) == p2_index);
        faces_adjacent_to_edge  = faces((faces_with_node1 & faces_with_node2), :);

        adj_face{i} = faces_adjacent_to_edge(~ismember(faces_adjacent_to_edge,[p1_index p2_index]));
end

% Only keep edges with 2 adjacent faces
mask = cellfun(@(x) numel(x)==2, adj_face);
creaseEdges = edges(mask,:);
creaseIdx   = find(mask);

% Build adjacency matrix safely
adj_mat = cell2mat(cellfun(@(x) x(:).', adj_face(mask), 'UniformOutput', false));

adj_faces.creaseIdx = creaseIdx;
adj_faces.adj_mat = adj_mat;
adj_faces.mask = mask;

% fix normals
focal_point = [0, 0, 99999];

p1    = vert_u(creaseEdges(:,1), :);
p2    = vert_u(creaseEdges(:,2), :);
p3_1  = vert_u(adj_mat(:,1), :);
p3_2  = vert_u(adj_mat(:,2), :);

e = p2 - p1;
e = e ./ vecnorm(e,2,2);

n1 = cross(p2 - p1, p3_1 - p1, 2);
n2 = cross(p3_2 - p1, p2 - p1, 2);

n1 = n1 ./ vecnorm(n1,2,2);
n2 = n2 ./ vecnorm(n2,2,2);

ref = focal_point-p1;
creaseEdges(dot(ref, n1, 2) < 0, :) = flip(creaseEdges(dot(ref, n1, 2) < 0, :), 2);
adj_faces.creaseEdges = creaseEdges;

end
