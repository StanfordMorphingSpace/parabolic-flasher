function angles = foldedCreaseAngles_fast(nodes_f, nodes, edges, adj_faces)
    nEdges = size(edges, 1);    
    angles = nan(nEdges, 1);
    
    % Only keep edges with 2 adjacent faces
    creaseEdges = adj_faces.creaseEdges;
    adj_mat = adj_faces.adj_mat;
    mask = adj_faces.mask;

    p1    = nodes_f(creaseEdges(:,1), :);
    p2    = nodes_f(creaseEdges(:,2), :);
    p3_1  = nodes_f(adj_mat(:,1), :);
    p3_2  = nodes_f(adj_mat(:,2), :);

    [temp_angles, ~, ~] = dihedralAngle(p1, p2, p3_1, p3_2);
    angles(mask) = real(temp_angles);
    
    % p1_u  = nodes(creaseEdges(:,1), :);
    % p2_u  = nodes(creaseEdges(:,2), :);
    % p3_1_u= nodes(adj_mat(:,1), :);
    % p3_2_u= nodes(adj_mat(:,2), :);
    % 
    % % Normals undeformed
    % n1_u = cross(p1_u-p3_1_u, p2_u-p3_1_u); n1_u = n1_u ./ vecnorm(n1_u,2,2);
    % n2_u = cross(p1_u-p3_2_u, p2_u-p3_2_u); n2_u = n2_u ./ vecnorm(n2_u,2,2);
    % 
    % % Normals deformed
    % n1 = cross(p1-p3_1, p2-p3_1); n1 = n1 ./ vecnorm(n1,2,2);
    % n2 = cross(p3_2-p1, p3_2-p2); n2 = n2 ./ vecnorm(n2,2,2);
    % 
    % % Flip orientation if needed
    % % n1(any(n1_u < 0, 2), :) = -n1(any(n1_u < 0, 2), :);
    % % n2(any(n2_u < 0, 2), :) = -n2(any(n2_u < 0, 2), :);
    % 
    % % Folding angles
    % temp_angles = angleBetweenVectors3d(n1,n2);
    % % flipMask = dot(n1,p3_2-p3_1,2)<0;
    % % temp_angles(flipMask) = -temp_angles(flipMask);
    % 
    % angles(mask) = real(temp_angles);
end