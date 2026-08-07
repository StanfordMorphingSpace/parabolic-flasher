function b = getb_f(nodes_f, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    h = geo.h;
    c = geo.c;
    n = geo.n;
    A = geo.A;
    cone_idx = varargin{1}{1}{2};
    stat_idxs = varargin{1}{1}{3};
    rib_n = varargin{1}{1}{4};
    
    surf_const_idxs = stat_idxs;
    if rib_n>0
       surf_const_idxs((end-rib_n+2):end) = [];
    end
    surf_const_idxs = surf_const_idxs(isnan(labels(surf_const_idxs)));
    
    % exclude inner nodes from spiral
    second_major_valley_indexes = find(labels < -1)';
    major_valley_indexes = arrayfun(@(x) find(labels == ...
        -labels(x, 1), 1), second_major_valley_indexes);
    
    num_surface_constrs = length(surf_const_idxs);
    num_spiral_constraints = length(major_valley_indexes);
    num_rotational_constrs = numel(second_major_valley_indexes)+1;
    b = zeros(num_rotational_constrs+num_spiral_constraints+num_surface_constrs+3,1);
    
    % g_spiral
    fv = @(x, y) mod(beta*(sqrt(x.^2+y.^2)-A)./(2*h), 2*pi) - mod(atan2(y, x), 2*pi); % valley
    
    spiral_constr_indexes = 1:num_spiral_constraints;
    x_coord_indexes = major_valley_indexes*3-2;
    y_coord_indexes = x_coord_indexes+1;
    b(spiral_constr_indexes) = fv(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));
    
    % add back inner nodes
    second_major_valley_indexes = find(labels < 0)';
    major_valley_indexes = arrayfun(@(x) find(labels == ...
        -labels(x, 1), 1), second_major_valley_indexes);
    
    % g_rot
    b(num_spiral_constraints + (1:3*num_rotational_constrs)) = getRotationalConstraints(major_valley_indexes, second_major_valley_indexes, nodes_f, beta);
    
    % g_inner (Applied as g_surface)
    surface_constraint_indexes = num_spiral_constraints + 3*num_rotational_constrs  + (1:(num_surface_constrs));
    b(surface_constraint_indexes)=  getSurfaceConstraints(surf_const_idxs, nodes_f, c, 0);

    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    z_coord_indexes = x_coord_indexes+2;
    b(end-2) = nodes_f(x_coord_indexes) - A;
    b(end-1) = nodes_f(y_coord_indexes);
    b(end) = nodes_f(z_coord_indexes) - c*A.^2;
end