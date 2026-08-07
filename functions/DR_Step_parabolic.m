function [p_u, p_f, v_u, v_f, a_u, a_f, E_crease, E_axial, E_v_u, E_v_f] = DR_Step_parabolic(p_u, p_f, v_u, v_f, ref, labels, edges, adj_faces, dt, geo, mass, i, stat_idxs, outer_idx, outer_R, cone_idx, n, rib_n, rib_d)
    beta    = 2*pi()/geo.N;

    angles_u = foldedCreaseAngles_fast(p_u, ref, edges, adj_faces);
    angles_f = foldedCreaseAngles_fast(p_f, ref, edges, adj_faces);
    lengths_u     = getEdgeLengths(p_u, edges);

    lengths_f   = getEdgeLengths(p_f, edges);

    % deployed
    [F_axial_u, ~, F_crease_u, ~, F_damping_u] = getDRForces_fast('u', p_u, ref, edges, adj_faces, geo.k_axial, geo.k_fold, lengths_f, angles_f, geo.gamma, v_u, mass);

    a_u = (F_axial_u + F_crease_u + F_damping_u)./mass;
    % a_u = (F_axial_u + F_crease_u)./mass;
    
    % folded
    [F_axial_f, E_axial, F_crease_f, E_crease, F_damping_f] = getDRForces_fast('f', p_f, ref, edges, adj_faces, geo.k_axial, geo.k_fold, lengths_u, angles_u, geo.gamma, v_f, mass);

    % calculate mass based on stiffnesses
    k_axial = geo.k_axial./lengths_u;
    k_crease = geo.k_fold.*lengths_u;

    a_f = (F_axial_f + F_damping_f + F_crease_f)./mass;
    % a_f = (F_axial_f + F_crease_f)./mass;

    for j = 1:length(stat_idxs) % exclude inner points so they stay fixed
        if j <= 2 || (rib_d > 0 && j > length(stat_idxs) - rib_n + 1)
            stat_a = [0,0,0];            
        else
            stat_a = (a_f(stat_idxs(j), :) + a_u(stat_idxs(j), :))./2;
        end
        a_u(stat_idxs(j), :) = stat_a; %stat_a; 
        a_f(stat_idxs(j), :) = stat_a; %stat_a;
    end
    
    % integrate with constraints
    [p_u, v_u, E_v_u] = makeStep(p_u, v_u, a_u, @getJacobian_u, @getb_u, labels, beta, i, mass, dt, {geo, outer_idx, stat_idxs, rib_d, rib_n});
    [p_f, v_f, E_v_f] = makeStep(p_f, v_f, a_f, @getJacobian_f, @getb_f, labels, beta, i, mass, dt, {geo, cone_idx, stat_idxs, rib_n});
    
end