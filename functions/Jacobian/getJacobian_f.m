function J = getJacobian_f(nodes_f, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    h = geo.h;
    c = geo.c;
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
    % g_spiral
    dfdx = @(x, y) y./(x.^2+y.^2) + beta*x./(2*h*sqrt(x.^2 + y.^2));
    dfdy = @(x, y) -x./(x.^2+y.^2) + beta*y./(2*h*sqrt(x.^2 + y.^2));
    
    J = sparse(num_rotational_constrs+num_spiral_constraints+num_surface_constrs+3, length(nodes_f)); % the +3 is for the fixed inner node
    
    spiral_constr_indexes = 1:num_spiral_constraints;
    x_coord_indexes = major_valley_indexes*3-2;
    y_coord_indexes = x_coord_indexes+1;
    J(sub2ind(size(J), spiral_constr_indexes, x_coord_indexes)) = dfdx(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));
    J(sub2ind(size(J), spiral_constr_indexes, y_coord_indexes)) = dfdy(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));
    
    % add back inner nodes
    second_major_valley_indexes = find(labels < 0)';
    major_valley_indexes = arrayfun(@(x) find(labels == ...
        -labels(x, 1), 1), second_major_valley_indexes);
    
    % g_rot
    rotational_constr_indexes = num_spiral_constraints + (1:3*num_rotational_constrs);
    Jrot = getRotationalConstraintJacobian(major_valley_indexes, ...
        second_major_valley_indexes, nodes_f, beta);
    J(rotational_constr_indexes, :) = Jrot;
    
    % g_inner (applied as a surface constraint)
    surface_constraint_indexes = num_spiral_constraints + 3*num_rotational_constrs  + (1:(num_surface_constrs));
    J(surface_constraint_indexes, :) = getSurfaceConstraintJacobian(surf_const_idxs, nodes_f, c);
    
    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    z_coord_indexes = x_coord_indexes+2;
    J(end-2, x_coord_indexes) = 1;
    J(end-1, y_coord_indexes) = 1;
    J(end, z_coord_indexes) = 1;
end