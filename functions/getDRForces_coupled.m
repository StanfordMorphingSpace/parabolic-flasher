function [F_axial_u, E_axial_u, F_crease_u, E_crease_u, F_damping_u, ...
          F_axial_f, E_axial_f, F_crease_f, E_crease_f, F_damping_f] = ...
    getDRForces_coupled(p_u, p_f, edges, adj_faces, EA, k_fold, gamma, v_u, v_f, mass)
    % Computes DR forces for the unfolded (u) and folded (f) networks together.
    % Each network's own edge lengths and crease dihedral angles are the other
    % network's rest lengths/angles.
    
    p_j_index = edges(:,1);
    p_k_index = edges(:,2);
    
    creaseEdges = adj_faces.creaseEdges;
    creaseIdx   = adj_faces.creaseIdx;
    adj_mat     = adj_faces.adj_mat;
    
    [edge_len_u, dir_u] = edgeLengthsAndDirs(p_u, p_j_index, p_k_index);
    [edge_len_f, dir_f] = edgeLengthsAndDirs(p_f, p_j_index, p_k_index);
    
    angle_u = creaseDihedralAngles(p_u, creaseEdges, adj_mat);
    angle_f = creaseDihedralAngles(p_f, creaseEdges, adj_mat);
    
    nEdges = size(edges, 1);
    angles_u_full = nan(nEdges, 1);
    angles_u_full(creaseIdx) = real(angle_u);
    angles_f_full = nan(nEdges, 1);
    angles_f_full(creaseIdx) = real(angle_f);
    
    [F_axial_u, E_axial_u, F_crease_u, E_crease_u, F_damping_u] = stateForces( ...
        'u', p_u, edge_len_u, dir_u, angle_u, edges, adj_faces, EA, k_fold, edge_len_f, angles_f_full, gamma, v_u, mass);
    
    [F_axial_f, E_axial_f, F_crease_f, E_crease_f, F_damping_f] = stateForces( ...
        'f', p_f, edge_len_f, dir_f, angle_f, edges, adj_faces, EA, k_fold, edge_len_u, angles_u_full, gamma, v_f, mass);
    
    end
    
function [edge_len, dir12] = edgeLengthsAndDirs(nodes_f, p_j_index, p_k_index)
    edge_vec = nodes_f(p_k_index, :) - nodes_f(p_j_index, :);
    edge_len = sqrt(sum(edge_vec.^2, 2));
    dir12 = edge_vec ./ edge_len;
end
    
function angle = creaseDihedralAngles(nodes_f, creaseEdges, adj_mat)
    p_j = nodes_f(creaseEdges(:,1), :);
    p_k = nodes_f(creaseEdges(:,2), :);
    p_i = nodes_f(adj_mat(:,1), :);
    p_l = nodes_f(adj_mat(:,2), :);
    angle = dihedralAngle(p_j, p_k, p_i, p_l);
end
    
function [F_axial, E_axial, F_crease, E_crease, F_damping] = stateForces(state, nodes_f, edge_len, dir12, angle, edges, adj_faces, EA, k_fold, lengths, angles, gamma, v, mass)
    % Same force/energy model as getDRForces_fast, but edge_len/dir12 (the own
    % network's edge geometry) and angle (its own crease dihedral angles) are
    % passed in already computed rather than recomputed from nodes_f.
    nEdges = length(edges);
    nNodes = length(nodes_f);
    F_crease = zeros(nNodes, 3);
    E_crease = zeros(nEdges, 1);
    
    p_j_index = edges(:,1);
    p_k_index = edges(:,2);
    
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
    
    rij = p_i-p_j;
    rkj = p_k-p_j;
    rkl = p_k-p_l;
    
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
    
    E_crease = zeros(nEdges,1);
    if state == 'u'
        E_crease(creaseIdx) = 1./2.*edge_len(creaseIdx).*k_fold(creaseIdx) .* dif_angles.^2;
    else
        E_crease(creaseIdx) = 1./2.*lengths(creaseIdx).*k_fold(creaseIdx) .* dif_angles.^2;
    end
end
