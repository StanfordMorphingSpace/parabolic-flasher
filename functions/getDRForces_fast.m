function [F_axial, E_axial, F_crease, E_crease, F_damping] = getDRForces_fast(state, nodes_f, nodes, edges, adj_faces, EA, k_fold, lengths, angles, gamma, v, mass)
    % nodes only used to calculate correct norm direction (based on flat)
    nEdges = length(edges);
    nNodes = length(nodes_f);
    F_crease = zeros(nNodes, 3);
    E_crease = zeros(nEdges, 1);

    p_j_index = edges(:,1);
    p_k_index = edges(:,2);

    p_j = nodes_f(p_j_index, :);
    p_k = nodes_f(p_k_index, :);
    
    % stretching
    edge_vec = p_k - p_j;
    edge_len = sqrt(sum(edge_vec.^2,2));
    dir12 = edge_vec ./ edge_len; 

    if state == 'u'
       coeff = (EA./edge_len).*(edge_len - lengths) - (EA./2./edge_len.^2).*(edge_len - lengths).^2;
    else
       coeff = (EA./lengths).*(edge_len - lengths);
    end

    F_axial = sparse(p_j_index, ones(nEdges,1), coeff .* dir12(:, 1), nNodes, 1) ...
        - sparse(p_k_index, ones(nEdges,1), coeff .* dir12(:, 1), nNodes, 1);
    F_axial(:,2) = sparse(p_j_index, ones(nEdges,1), coeff .* dir12(:, 2), nNodes, 1) ...
                  - sparse(p_k_index, ones(nEdges,1), coeff .* dir12(:, 2), nNodes, 1);
    F_axial(:,3) = sparse(p_j_index, ones(nEdges,1), coeff .* dir12(:, 3), nNodes, 1) ...
                  - sparse(p_k_index, ones(nEdges,1), coeff .* dir12(:, 3), nNodes, 1);
    F_axial = full(F_axial);


    if state == 'u'
        E_axial = 0.5*(EA./edge_len).*(edge_len - lengths).^2;
    else
        E_axial = 0.5*(EA./lengths).*(edge_len - lengths).^2;
    end
    

    % damping
    if state == 'u'
        c1 = 2*gamma*sqrt(EA./edge_len .* mass(p_j_index));
        c2 = 2*gamma*sqrt(EA./edge_len .* mass(p_k_index));
    else
        c1 = 2*gamma*sqrt(EA./lengths .* mass(p_j_index));
        c2 = 2*gamma*sqrt(EA./lengths .* mass(p_k_index));
    end


    dv = v(p_k_index,:) - v(p_j_index,:);

    F_damping = sparse(p_j_index, ones(nEdges,1), c1 .* dv(:, 1), nNodes, 1) ...
                  - sparse(p_k_index, ones(nEdges,1), c2 .* dv(:, 1), nNodes, 1);
    F_damping(:,2) = sparse(p_j_index, ones(nEdges,1), c1 .* dv(:, 2), nNodes, 1) ...
                  - sparse(p_k_index, ones(nEdges,1), c2 .* dv(:, 2), nNodes, 1);
    F_damping(:,3) = sparse(p_j_index, ones(nEdges,1), c1 .* dv(:, 3), nNodes, 1) ...
                  - sparse(p_k_index, ones(nEdges,1), c2 .* dv(:, 3), nNodes, 1);
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

    p_j    = nodes_f(creaseEdges(:,1), :);
    p_k    = nodes_f(creaseEdges(:,2), :);
    p_i  = nodes_f(adj_mat(:,1), :);
    p_l  = nodes_f(adj_mat(:,2), :);

    [angle, m, n] = dihedralAngle(p_j, p_k, p_i, p_l);
    % angle = real(angle);

    dif_angles = angle - angles(creaseIdx);

    % stiffness
    if state == 'u'
        scale = -edge_len(creaseIdx).*k_fold(creaseIdx) .* dif_angles;
        dif_k_crease = zeros(nEdges,1);
        dif_k_crease(creaseIdx) = 1./2.*k_fold(creaseIdx) .* dif_angles.^2;
        eq_axial = sparse(p_j_index, ones(nEdges,1), dif_k_crease .* dir12(:, 1), nNodes, 1) ...
                 - sparse(p_k_index, ones(nEdges,1), dif_k_crease .* dir12(:, 1), nNodes, 1);
        eq_axial(:,2) = sparse(p_j_index, ones(nEdges,1), dif_k_crease .* dir12(:, 2), nNodes, 1) ...
                      - sparse(p_k_index, ones(nEdges,1), dif_k_crease .* dir12(:, 2), nNodes, 1);
        eq_axial(:,3) = sparse(p_j_index, ones(nEdges,1), dif_k_crease .* dir12(:, 3), nNodes, 1) ...
                      - sparse(p_k_index, ones(nEdges,1), dif_k_crease .* dir12(:, 3), nNodes, 1);
        eq_axial = full(eq_axial);
    else
        scale = -lengths(creaseIdx).*k_fold(creaseIdx) .* dif_angles;
        eq_axial = zeros(size(F_crease));
    end
    
    
    % Forces
    %dif_angles = mod((angle - angles(creaseIdx)) + pi, 2*pi) - pi;

    rij = p_i-p_j;
    rkj = p_k-p_j;
    rkl = p_k-p_l;

    % Raw normals
    m = cross(rij, rkj, 2);
    n = cross(rkj, rkl, 2);

    norm_kj = vecnorm(rkj, 2, 2).^2;

    dthdp_i = vecnorm(rkj, 2, 2)./vecnorm(m, 2, 2).^2 .* m;
    dthdp_l = -vecnorm(rkj, 2, 2)./vecnorm(n, 2, 2).^2 .* n;
    dthdp_j = (sum(rij.*rkj,2)./norm_kj - 1).*dthdp_i - sum(rkl.*rkj,2)./norm_kj.*dthdp_l;
    dthdp_k = (sum(rkl.*rkj,2)./norm_kj - 1).*dthdp_l - sum(rij.*rkj,2)./norm_kj.*dthdp_i;


    F_crease = sparse(adj_mat(:,1), ones(nCrease,1), scale.*dthdp_i(:, 1), nNodes, 1) ...
                  + sparse(adj_mat(:,2), ones(nCrease,1), scale.*dthdp_l(:, 1), nNodes, 1) ...
                  + sparse(creaseEdges(:,1), ones(nCrease,1), scale.*dthdp_j(:, 1), nNodes, 1) ...
                  + sparse(creaseEdges(:,2), ones(nCrease,1), scale.*dthdp_k(:, 1), nNodes, 1);
    F_crease(:,2) = sparse(adj_mat(:,1), ones(nCrease,1), scale.*dthdp_i(:, 2), nNodes, 1) ...
                  + sparse(adj_mat(:,2), ones(nCrease,1), scale.*dthdp_l(:, 2), nNodes, 1) ...
                  + sparse(creaseEdges(:,1), ones(nCrease,1), scale.*dthdp_j(:, 2), nNodes, 1) ...
                  + sparse(creaseEdges(:,2), ones(nCrease,1), scale.*dthdp_k(:, 2), nNodes, 1);
    F_crease(:,3) = sparse(adj_mat(:,1), ones(nCrease,1), scale.*dthdp_i(:, 3), nNodes, 1) ...
                  + sparse(adj_mat(:,2), ones(nCrease,1), scale.*dthdp_l(:, 3), nNodes, 1) ...
                  + sparse(creaseEdges(:,1), ones(nCrease,1), scale.*dthdp_j(:, 3), nNodes, 1) ...
                  + sparse(creaseEdges(:,2), ones(nCrease,1), scale.*dthdp_k(:, 3), nNodes, 1);
    F_crease = full(F_crease);

    F_crease = F_crease + eq_axial;
    
    % Energy (zeros elsewhere)
    E_crease = zeros(nEdges,1);
    if state == 'u'
        E_crease(creaseIdx) = 1./2.*edge_len(creaseIdx).*k_fold(creaseIdx) .* dif_angles.^2;
    else
        E_crease(creaseIdx) = 1./2.*lengths(creaseIdx).*k_fold(creaseIdx) .* dif_angles.^2;
    end

    % %% old version
    % 
    % p_j_u  = nodes(creaseEdges(:,1), :);
    % p_k_u  = nodes(creaseEdges(:,2), :);
    % p_i_u= nodes(adj_mat(:,1), :);
    % p_l_u= nodes(adj_mat(:,2), :);
    % % 
    % % Heights
    % h1 = vecnorm(cross(p_j-p_i, p_k-p_i),2,2)./vecnorm(p_k-p_j,2,2);
    % h2 = vecnorm(cross(p_j-p_l, p_k-p_l),2,2)./vecnorm(p_k-p_j,2,2);
    % 
    % % Normals undeformed
    % n1_u = cross(p_j_u-p_i_u, p_k_u-p_i_u); n1_u = n1_u ./ vecnorm(n1_u,2,2);
    % n2_u = cross(p_j_u-p_l_u, p_k_u-p_l_u); n2_u = n2_u ./ vecnorm(n2_u,2,2);
    % 
    % % Normals deformed
    % n1 = cross(p_k - p_j, p_i - p_j, 2);
    % n2 = cross(p_l - p_j, p_k - p_j, 2);
    % 
    % % Normalize
    % n1 = n1 ./ vecnorm(n1,2,2);
    % n2 = n2 ./ vecnorm(n2,2,2);
    % 
    % % Interior angles
    % a431 = angleBetweenVectors3d(p_j-p_k,  p_i-p_k);
    % a314 = angleBetweenVectors3d(p_i-p_j, p_k-p_j);
    % a423 = angleBetweenVectors3d(p_l-p_k, p_j-p_k);
    % a342 = angleBetweenVectors3d(p_k-p_j,  p_l-p_j);
    % 
    % dthdp_j = -cot(a431)./(cot(a314)+cot(a431)).*n1./h1 + ...
    %          -cot(a423)./(cot(a342)+cot(a423)).*n2./h2;
    % dthdp_k = -cot(a314)./(cot(a314)+cot(a431)).*n1./h1 + ...
    %          -cot(a342)./(cot(a342)+cot(a423)).*n2./h2;
    % 
    % 
    % 
    % F_crease2 = sparse(adj_mat(:,1), ones(nCrease,1), -scale.*n1(:, 1)./h1, nNodes, 1) ...
    %               + sparse(adj_mat(:,2), ones(nCrease,1), -scale.*n2(:, 1)./h2, nNodes, 1) ...
    %               + sparse(creaseEdges(:,1), ones(nCrease,1), -scale.*dthdp_j(:, 1), nNodes, 1) ...
    %               + sparse(creaseEdges(:,2), ones(nCrease,1), -scale.*dthdp_k(:, 1), nNodes, 1);
    % F_crease2(:,2) = sparse(adj_mat(:,1), ones(nCrease,1), -scale.*n1(:, 2)./h1, nNodes, 1) ...
    %               + sparse(adj_mat(:,2), ones(nCrease,1), -scale.*n2(:, 2)./h2, nNodes, 1) ...
    %               + sparse(creaseEdges(:,1), ones(nCrease,1), -scale.*dthdp_j(:, 2), nNodes, 1) ...
    %               + sparse(creaseEdges(:,2), ones(nCrease,1), -scale.*dthdp_k(:, 2), nNodes, 1);
    % F_crease2(:,3) = sparse(adj_mat(:,1), ones(nCrease,1), -scale.*n1(:, 3)./h1, nNodes, 1) ...
    %               + sparse(adj_mat(:,2), ones(nCrease,1), -scale.*n2(:, 3)./h2, nNodes, 1) ...
    %               + sparse(creaseEdges(:,1), ones(nCrease,1), -scale.*dthdp_j(:, 3), nNodes, 1) ...
    %               + sparse(creaseEdges(:,2), ones(nCrease,1), -scale.*dthdp_k(:, 3), nNodes, 1);
    % F_crease2 = full(F_crease2);

    %F_crease = F_crease2 + eq_axial;
end