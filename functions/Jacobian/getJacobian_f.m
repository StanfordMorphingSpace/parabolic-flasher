function J = getJacobian_f(nodes_f, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    h = geo.h;
    surf_func_prime = geo.surf_func_prime;

    stat_idxs = labels.stat_idxs;
    valley_idx = labels.valley_idx;
    valley_rot_idx = labels.valley_rot_idx;

    % the brim continues the valley/rotated-valley lines beyond R onto the
    % same folded spiral, so it's just more pairs in the same constraints
    has_brim = isfield(labels, 'brim_idx');
    if has_brim
        valley_idx = [valley_idx, labels.brim_valley_idx];
        valley_rot_idx = [valley_rot_idx, labels.brim_valley_rot_idx];
    end

    % the inner-edge nodes shared between the two gores are pinned to the
    % surface directly instead of via the spiral/rotational constraints;
    % excludes the rib/brim attachment points (not on this surface) and the fixed inner node
    exclude_from_surf = [valley_idx, valley_rot_idx];
    if isfield(labels, 'ribs_idx')
        exclude_from_surf = [exclude_from_surf, labels.ribs_idx];
    end
    if has_brim
        exclude_from_surf = [exclude_from_surf, labels.brim_idx];
    end
    surf_const_idxs = stat_idxs(~ismember(stat_idxs, exclude_from_surf));

    % exclude the fixed inner node from the spiral (constrained separately)
    spiral_major = valley_idx(2:(end-length(labels.brim_valley_idx)));

    num_surface_constrs = length(surf_const_idxs);
    num_spiral_constraints = length(spiral_major);
    num_rotational_constrs = length(valley_rot_idx);

    % g_spiral
    dfdx = @(x, y) y./(x.^2+y.^2) + beta*x./(2*h*sqrt(x.^2 + y.^2));
    dfdy = @(x, y) -x./(x.^2+y.^2) + beta*y./(2*h*sqrt(x.^2 + y.^2));

    J = sparse(num_rotational_constrs+num_spiral_constraints+num_surface_constrs+3, length(nodes_f)); % the +3 is for the fixed inner node

    spiral_constr_indexes = 1:num_spiral_constraints;
    x_coord_indexes = spiral_major*3-2;
    y_coord_indexes = x_coord_indexes+1;
    J(sub2ind(size(J), spiral_constr_indexes, x_coord_indexes)) = dfdx(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));
    J(sub2ind(size(J), spiral_constr_indexes, y_coord_indexes)) = dfdy(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));

    % g_rot
    rotational_constr_indexes = num_spiral_constraints + (1:3*num_rotational_constrs);
    Jrot = getRotationalConstraintJacobian(valley_idx, valley_rot_idx, nodes_f, beta);
    J(rotational_constr_indexes, :) = Jrot;

    % g_inner (applied as a surface constraint)
    surface_constraint_indexes = num_spiral_constraints + 3*num_rotational_constrs  + (1:(num_surface_constrs));
    J(surface_constraint_indexes, :) = getSurfaceConstraintJacobian(surf_const_idxs, nodes_f, surf_func_prime);

    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    z_coord_indexes = x_coord_indexes+2;
    J(end-2, x_coord_indexes) = 1;
    J(end-1, y_coord_indexes) = 1;
    J(end, z_coord_indexes) = 1;
end
