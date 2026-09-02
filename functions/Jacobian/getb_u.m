function b = getb_u(nodes_u, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    outer_idx = labels.outer_idx;
    stat_idxs = labels.stat_idxs;
    valley_idx = labels.valley_idx;
    valley_rot_idx = labels.valley_rot_idx;
    num_outer_nodes = length(outer_idx);
    surf_func = geo.surf_func;
    A = geo.A;
    outer_R = geo.R;
    rib_d = geo.rib_d;

    has_brim = isfield(labels, 'brim_idx');
    if has_brim
        brim_func = geo.brim_func;
        brim_outer_idx = labels.brim_outer_idx;
        valley_idx = [valley_idx, labels.brim_valley_idx];
        valley_rot_idx = [valley_rot_idx, labels.brim_valley_rot_idx];
    else
        brim_outer_idx = [];
    end
    num_brim_outer_nodes = length(brim_outer_idx);

    num_rot = length(valley_rot_idx);
    nNodes = length(nodes_u)/3;

    has_ribs = isfield(labels, 'ribs_idx');
    if has_ribs
        ribs_lower_idx = labels.ribs_lower_idx;
    else
        ribs_lower_idx = [];
    end

    exclude_from_surf = valley_rot_idx;
    if has_ribs
        exclude_from_surf = [exclude_from_surf, labels.ribs_idx];
    end
    if has_brim
        exclude_from_surf = [exclude_from_surf, labels.brim_idx];
        brim_surface_node_indexes = reshape(setdiff(labels.brim_idx, labels.brim_valley_rot_idx), 1, []);
    else
        brim_surface_node_indexes = [];
    end
    surface_node_indexes = reshape(setdiff(1:nNodes, exclude_from_surf), 1, []);
    num_brim_surf = length(brim_surface_node_indexes);

    b = zeros(num_outer_nodes + num_brim_outer_nodes + length(surface_node_indexes) + length(ribs_lower_idx) + num_brim_surf + num_rot*3 + 2, 1);

    % g_outer
    x_coord_indexes = outer_idx.*3-2;
    y_coord_indexes = x_coord_indexes+1;
    b(1:num_outer_nodes)= -outer_R^2 + (nodes_u(x_coord_indexes).^2 + nodes_u(y_coord_indexes).^2);

    % g_brim_outer
    if has_brim
        x_coord_indexes = brim_outer_idx.*3-2;
        y_coord_indexes = x_coord_indexes+1;
        b(num_outer_nodes+(1:num_brim_outer_nodes)) = -geo.brim_R^2 + (nodes_u(x_coord_indexes).^2 + nodes_u(y_coord_indexes).^2);
    end

    % g_rot
    b(num_outer_nodes+num_brim_outer_nodes+(1:3*num_rot)) = getRotationalConstraints(valley_idx, valley_rot_idx, nodes_u, beta);

    % g_surface - everyone except the rotated valley copies, brim nodes, and (all) rib nodes
    surface_constraint_indexes = num_rot * 3 + num_outer_nodes + num_brim_outer_nodes + (1:length(surface_node_indexes));
    b(surface_constraint_indexes) = getSurfaceConstraints(surface_node_indexes, nodes_u, surf_func, 0);

    % g_surface (rib) - the rib's bottom row is pinned to a surface offset by -rib_d
    if has_ribs
        rib_constraint_indexes = num_rot * 3 + num_outer_nodes + num_brim_outer_nodes + length(surface_node_indexes) + (1:length(ribs_lower_idx));
        b(rib_constraint_indexes) = getSurfaceConstraints(ribs_lower_idx, nodes_u, surf_func, -rib_d);
    end

    % g_brim_surface - brim nodes (except their rotated-valley copies) stay on brim_func
    if has_brim
        brim_surface_constraint_indexes = num_rot * 3 + num_outer_nodes + num_brim_outer_nodes + length(surface_node_indexes) + length(ribs_lower_idx) + (1:num_brim_surf);
        b(brim_surface_constraint_indexes) = getSurfaceConstraints(brim_surface_node_indexes, nodes_u, brim_func, 0);
    end

    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    b(end-1) = nodes_u(x_coord_indexes) - A;
    b(end) = nodes_u(y_coord_indexes);

end
