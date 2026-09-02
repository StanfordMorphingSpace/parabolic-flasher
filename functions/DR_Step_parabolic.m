function [p_u, p_f, v_u, v_f, a_u, a_f, E_crease, E_axial, E_v_u, E_v_f] = DR_Step_parabolic(p_u, p_f, v_u, v_f, labels, edges, adj_faces, dt, geo, mass, i)
    beta    = 2*pi()/geo.N;
    stat_idxs = labels.stat_idxs;

    [F_axial_u, ~, F_crease_u, ~, F_damping_u, ...
     F_axial_f, E_axial, F_crease_f, E_crease, F_damping_f] = ...
        getDRForces_coupled(p_u, p_f, edges, adj_faces, geo.k_axial, geo.k_fold, geo.gamma, v_u, v_f, mass);

    a_u = (F_axial_u + F_crease_u + F_damping_u)./mass;
    a_f = (F_axial_f + F_crease_f + F_damping_f)./mass;

    % rib attachment points (the trailing entries of stat_idxs, one per rib row)
    % are fully fixed, same as the two hub nodes; everything else is averaged
    has_ribs = isfield(labels, 'ribs_idx');
    if has_ribs
        num_rib_attach = sum(ismember(stat_idxs, labels.ribs_idx));
    else
        num_rib_attach = 0;
    end

    for j = 1:length(stat_idxs) % exclude inner points so they stay fixed
        if j <= 2 || (has_ribs && j > length(stat_idxs) - num_rib_attach)
            stat_a = [0,0,0];
        else
            stat_a = (a_f(stat_idxs(j), :) + a_u(stat_idxs(j), :))./2;
        end
        a_u(stat_idxs(j), :) = stat_a; %stat_a;
        a_f(stat_idxs(j), :) = stat_a; %stat_a;
    end

    % integrate with constraints
    [p_u, v_u, E_v_u] = makeStep(p_u, v_u, a_u, @getJacobian_u, @getb_u, labels, beta, i, mass, dt, {geo});
    [p_f, v_f, E_v_f] = makeStep(p_f, v_f, a_f, @getJacobian_f, @getb_f, labels, beta, i, mass, dt, {geo});

end