function [F_axial, E_axial, F_crease, E_crease, F_damping] = getDRForces_fast(nodes_f, nodes, edges, adj_faces, EA, k_fold, lengths, angles, gamma, v, mass)
    % nodes only used to calculate correct norm direction (based on flat)
    nEdges = length(edges);
    nNodes = length(nodes_f);
    F_crease = zeros(nNodes, 3);
    E_crease = zeros(nEdges, 1);

    p1_index = edges(:,1);
    p2_index = edges(:,2);

    p1 = nodes_f(p1_index, :);
    p2 = nodes_f(p2_index, :);
    
    % stretching
    edge_vec = p2 - p1;
    edge_len = sqrt(sum(edge_vec.^2,2));
    dir12 = edge_vec ./ edge_len; 

    coeff = (EA./lengths).*(edge_len - lengths);

    F_axial = sparse(p1_index, ones(nEdges,1), coeff .* dir12(:, 1), nNodes, 1) ...
        - sparse(p2_index, ones(nEdges,1), coeff .* dir12(:, 1), nNodes, 1);
    F_axial(:,2) = sparse(p1_index, ones(nEdges,1), coeff .* dir12(:, 2), nNodes, 1) ...
                  - sparse(p2_index, ones(nEdges,1), coeff .* dir12(:, 2), nNodes, 1);
    F_axial(:,3) = sparse(p1_index, ones(nEdges,1), coeff .* dir12(:, 3), nNodes, 1) ...
                  - sparse(p2_index, ones(nEdges,1), coeff .* dir12(:, 3), nNodes, 1);
    F_axial = full(F_axial);

    E_axial = 0.5*(EA./lengths).*(edge_len - lengths).^2;

    % damping
    c1 = 2*gamma*sqrt(EA./lengths .* mass(p1_index));
    c2 = 2*gamma*sqrt(EA./lengths .* mass(p2_index));

    dv = v(p2_index,:) - v(p1_index,:);

    F_damping = sparse(p1_index, ones(nEdges,1), c1 .* dv(:, 1), nNodes, 1) ...
                  - sparse(p2_index, ones(nEdges,1), c2 .* dv(:, 1), nNodes, 1);
    F_damping(:,2) = sparse(p1_index, ones(nEdges,1), c1 .* dv(:, 2), nNodes, 1) ...
                  - sparse(p2_index, ones(nEdges,1), c2 .* dv(:, 2), nNodes, 1);
    F_damping(:,3) = sparse(p1_index, ones(nEdges,1), c1 .* dv(:, 3), nNodes, 1) ...
                  - sparse(p2_index, ones(nEdges,1), c2 .* dv(:, 3), nNodes, 1);
    F_damping = full(F_damping);

    % bending
    if all(angles == 0)
        F_crease = zeros(nNodes,3);
        E_crease = zeros(nEdges,1);
        return;
    end

    creaseEdges = adj_faces.creaseEdges;
    creaseIdx   = adj_faces.creaseIdx;
    adj_mat = adj_faces.adj_mat;

    nCrease = length(adj_mat);

    p1    = nodes_f(creaseEdges(:,1), :);
    p2    = nodes_f(creaseEdges(:,2), :);
    p3_1  = nodes_f(adj_mat(:,1), :);
    p3_2  = nodes_f(adj_mat(:,2), :);

    [angle, n1, n2] = dihedralAngle(p1, p2, p3_1, p3_2);
    angle = real(angle);

    
    p1_u  = nodes(creaseEdges(:,1), :);
    p2_u  = nodes(creaseEdges(:,2), :);
    p3_1_u= nodes(adj_mat(:,1), :);
    p3_2_u= nodes(adj_mat(:,2), :);
    
    % Heights
    h1 = vecnorm(cross(p1-p3_1, p2-p3_1),2,2)./vecnorm(p2-p1,2,2);
    h2 = vecnorm(cross(p1-p3_2, p2-p3_2),2,2)./vecnorm(p2-p1,2,2);
    
    % Normals undeformed
    % n1_u = cross(p1_u-p3_1_u, p2_u-p3_1_u); n1_u = n1_u ./ vecnorm(n1_u,2,2);
    % n2_u = cross(p1_u-p3_2_u, p2_u-p3_2_u); n2_u = n2_u ./ vecnorm(n2_u,2,2);
    
    % Normals deformed
    % n1 = cross(p1-p3_1, p2-p3_1); n1 = n1 ./ vecnorm(n1,2,2);
    % n2 = cross(p3_2-p1, p3_2-p2); n2 = n2 ./ vecnorm(n2,2,2);
    
    % Flip orientation if needed
    % n1(any(n1_u < 0, 2), :) = -n1(any(n1_u < 0, 2), :);
    % n2(any(n2_u < 0, 2), :) = -n2(any(n2_u < 0, 2), :);
    
    % Folding angles
    % angle = angleBetweenVectors3d(n1,n2);
    % angle = real(angle);
    % flipMask = dot(n1,p3_2-p3_1,2)<0;
    % angle(flipMask) = -angle(flipMask);
    
    % Interior angles
    a431 = angleBetweenVectors3d(p1-p2,  p3_1-p2);
    a314 = angleBetweenVectors3d(p3_1-p1, p2-p1);
    a423 = angleBetweenVectors3d(p3_2-p2, p1-p2);
    a342 = angleBetweenVectors3d(p2-p1,  p3_2-p1);
        
    dthdp1 = -cot(a431)./(cot(a314)+cot(a431)).*n1./h1 + ...
             -cot(a423)./(cot(a342)+cot(a423)).*n2./h2;
    dthdp2 = -cot(a314)./(cot(a314)+cot(a431)).*n1./h1 + ...
             -cot(a342)./(cot(a342)+cot(a423)).*n2./h2;
    
    % stiffness
    k_crease = lengths(creaseIdx).*k_fold(creaseIdx);
    
    % Forces
    %dif_angles = mod((angle - angles(creaseIdx)) + pi, 2*pi) - pi;
    dif_angles = angle - angles(creaseIdx);
    scale = k_crease .* dif_angles;

    F_crease = sparse(adj_mat(:,1), ones(nCrease,1), -scale.*n1(:, 1)./h1, nNodes, 1) ...
                  + sparse(adj_mat(:,2), ones(nCrease,1), -scale.*n2(:, 1)./h2, nNodes, 1) ...
                  + sparse(creaseEdges(:,1), ones(nCrease,1), -scale.*dthdp1(:, 1), nNodes, 1) ...
                  + sparse(creaseEdges(:,2), ones(nCrease,1), -scale.*dthdp2(:, 1), nNodes, 1);
    F_crease(:,2) = sparse(adj_mat(:,1), ones(nCrease,1), -scale.*n1(:, 2)./h1, nNodes, 1) ...
                  + sparse(adj_mat(:,2), ones(nCrease,1), -scale.*n2(:, 2)./h2, nNodes, 1) ...
                  + sparse(creaseEdges(:,1), ones(nCrease,1), -scale.*dthdp1(:, 2), nNodes, 1) ...
                  + sparse(creaseEdges(:,2), ones(nCrease,1), -scale.*dthdp2(:, 2), nNodes, 1);
    F_crease(:,3) = sparse(adj_mat(:,1), ones(nCrease,1), -scale.*n1(:, 3)./h1, nNodes, 1) ...
                  + sparse(adj_mat(:,2), ones(nCrease,1), -scale.*n2(:, 3)./h2, nNodes, 1) ...
                  + sparse(creaseEdges(:,1), ones(nCrease,1), -scale.*dthdp1(:, 3), nNodes, 1) ...
                  + sparse(creaseEdges(:,2), ones(nCrease,1), -scale.*dthdp2(:, 3), nNodes, 1);
    F_crease = full(F_crease);
    
    % Energy (zeros elsewhere)
    E_crease = zeros(nEdges,1);
    E_crease(creaseIdx) = 0.5*k_crease .* (angle - angles(creaseIdx)).^2;
end