function J = getJacobian_u(nodes_u, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    surf_func_prime = geo.surf_func_prime;

    outer_idx = labels.outer_idx;
    stat_idxs = labels.stat_idxs;
    valley_idx = labels.valley_idx;
    valley_rot_idx = labels.valley_rot_idx;

    % the brim continues the valley/rotated-valley lines beyond R, so its
    % rotational-periodicity constraint is just more pairs of the same kind
    has_brim = isfield(labels, 'brim_idx');
    if has_brim
        surf_func_prime_brim = geo.brim_func_prime;
        brim_outer_idx = labels.brim_outer_idx;
        valley_idx = [valley_idx, labels.brim_valley_idx];
        valley_rot_idx = [valley_rot_idx, labels.brim_valley_rot_idx];
    else
        brim_outer_idx = [];
    end

    num_outer_nodes = length(outer_idx);
    num_brim_outer_nodes = length(brim_outer_idx);
    num_rot = length(valley_rot_idx);
    nNodes = length(nodes_u)/3;

    % g_surface applies to everyone except the rotated valley copies, brim
    % nodes (which use brim_func instead, in their own block below), and
    % partial-depth rib rows (only the bottom row is pinned to a surface)
    exclude_from_surf = valley_rot_idx;
    if isfield(labels, 'ribs_idx')
        exclude_from_surf = [exclude_from_surf, setdiff(labels.ribs_idx, labels.ribs_lower_idx)];
    end
    if has_brim
        exclude_from_surf = [exclude_from_surf, labels.brim_idx];
        brim_surface_node_indexes = reshape(setdiff(labels.brim_idx, labels.brim_valley_rot_idx), 1, []);
    else
        brim_surface_node_indexes = [];
    end
    surface_node_indexes = reshape(setdiff(1:nNodes, exclude_from_surf), 1, []);
    num_brim_surf = length(brim_surface_node_indexes);

    J = sparse(num_outer_nodes + num_brim_outer_nodes + length(surface_node_indexes) + num_brim_surf + num_rot*3 + 2, length(nodes_u)); % the +2 is from fixed inner node

    % g_outer
    x_coord_indexes = outer_idx*3-2;
    y_coord_indexes = outer_idx*3-1;
    J(sub2ind(size(J),1:num_outer_nodes, x_coord_indexes)) = 2*nodes_u(x_coord_indexes); %2*dfdx
    J(sub2ind(size(J),1:num_outer_nodes, y_coord_indexes)) = 2*nodes_u(y_coord_indexes); %2*dfdy

    % g_brim_outer - same form as g_outer, for the brim's own outer radius
    if has_brim
        x_coord_indexes = brim_outer_idx*3-2;
        y_coord_indexes = brim_outer_idx*3-1;
        rows = num_outer_nodes + (1:num_brim_outer_nodes);
        J(sub2ind(size(J), rows, x_coord_indexes)) = 2*nodes_u(x_coord_indexes);
        J(sub2ind(size(J), rows, y_coord_indexes)) = 2*nodes_u(y_coord_indexes);
    end

    % g_rot
    rotational_constr_indexes = num_outer_nodes + num_brim_outer_nodes + (1:3*num_rot);
    Jrot = getRotationalConstraintJacobian(valley_idx, valley_rot_idx, nodes_u, beta);
    J(rotational_constr_indexes, :) = Jrot;

    % g_surface - apply to all but the rotated valley copies, brim nodes, and mid-depth rib rows
    surface_constraint_indexes = num_rot*3 + num_outer_nodes + num_brim_outer_nodes + (1:length(surface_node_indexes));
    J(surface_constraint_indexes, :) = getSurfaceConstraintJacobian(surface_node_indexes, nodes_u, surf_func_prime);

    % g_brim_surface - brim nodes (except their rotated-valley copies) stay on brim_func
    if has_brim
        brim_surface_constraint_indexes = num_rot*3 + num_outer_nodes + num_brim_outer_nodes + length(surface_node_indexes) + (1:num_brim_surf);
        J(brim_surface_constraint_indexes, :) = getSurfaceConstraintJacobian(brim_surface_node_indexes, nodes_u, surf_func_prime_brim);
    end

    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    J(end-1, x_coord_indexes) = 1;
    J(end, y_coord_indexes) = 1;
end
