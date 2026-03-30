function angles = foldedCreaseAngles_fast(nodes_f, nodes, edges, adj_faces)
    nEdges = size(edges, 1);    
    angles = nan(nEdges, 1);
    
    % Only keep edges with 2 adjacent faces
    creaseEdges = adj_faces.creaseEdges;
    adj_mat = adj_faces.adj_mat;
    mask = adj_faces.mask;

    p_j    = nodes_f(creaseEdges(:,1), :);
    p_k    = nodes_f(creaseEdges(:,2), :);
    p_i  = nodes_f(adj_mat(:,1), :);
    p_l  = nodes_f(adj_mat(:,2), :);

    [temp_angles, ~, ~] = dihedralAngle(p_j, p_k, p_i, p_l);
    angles(mask) = real(temp_angles);
    
    % p_j_u  = nodes(creaseEdges(:,1), :);
    % p_k_u  = nodes(creaseEdges(:,2), :);
    % p_i_u= nodes(adj_mat(:,1), :);
    % p_l_u= nodes(adj_mat(:,2), :);
    % 
    % % Normals undeformed
    % n1_u = cross(p_j_u-p_i_u, p_k_u-p_i_u); n1_u = n1_u ./ vecnorm(n1_u,2,2);
    % n2_u = cross(p_j_u-p_l_u, p_k_u-p_l_u); n2_u = n2_u ./ vecnorm(n2_u,2,2);
    % 
    % % Normals deformed
    % n1 = cross(p_j-p_i, p_k-p_i); n1 = n1 ./ vecnorm(n1,2,2);
    % n2 = cross(p_l-p_j, p_l-p_k); n2 = n2 ./ vecnorm(n2,2,2);
    % 
    % % Flip orientation if needed
    % % n1(any(n1_u < 0, 2), :) = -n1(any(n1_u < 0, 2), :);
    % % n2(any(n2_u < 0, 2), :) = -n2(any(n2_u < 0, 2), :);
    % 
    % % Folding angles
    % temp_angles = angleBetweenVectors3d(n1,n2);
    % % flipMask = dot(n1,p_l-p_i,2)<0;
    % % temp_angles(flipMask) = -temp_angles(flipMask);
    % 
    % angles(mask) = real(temp_angles);
end