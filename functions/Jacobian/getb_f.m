function b = getb_f(nodes_f, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    h = geo.h;
    surf_func = geo.surf_func;
    A = geo.A;

    stat_idxs = labels.stat_idxs;
    valley_idx = labels.valley_idx;
    valley_rot_idx = labels.valley_rot_idx;

    has_brim = isfield(labels, 'brim_idx');
    if has_brim
        valley_idx = [valley_idx, labels.brim_valley_idx];
        valley_rot_idx = [valley_rot_idx, labels.brim_valley_rot_idx];
    end

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
    b = zeros(num_rotational_constrs+num_spiral_constraints+num_surface_constrs+3,1);

    % g_spiral
    fv = @(x, y) mod(beta*(sqrt(x.^2+y.^2)-A)./(2*h), 2*pi) - mod(atan2(y, x), 2*pi); % valley

    spiral_constr_indexes = 1:num_spiral_constraints;
    x_coord_indexes = spiral_major*3-2;
    y_coord_indexes = x_coord_indexes+1;
    b(spiral_constr_indexes) = fv(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));

    % g_rot
    b(num_spiral_constraints + (1:3*num_rotational_constrs)) = getRotationalConstraints(valley_idx, valley_rot_idx, nodes_f, beta);

    % g_inner (Applied as g_surface)
    surface_constraint_indexes = num_spiral_constraints + 3*num_rotational_constrs  + (1:(num_surface_constrs));
    b(surface_constraint_indexes)=  getSurfaceConstraints(surf_const_idxs, nodes_f, surf_func, 0);

    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    z_coord_indexes = x_coord_indexes+2;
    b(end-2) = nodes_f(x_coord_indexes) - A;
    b(end-1) = nodes_f(y_coord_indexes);
    b(end) = nodes_f(z_coord_indexes) - surf_func(A);
end
