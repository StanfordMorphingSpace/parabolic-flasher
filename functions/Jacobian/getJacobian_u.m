function J = getJacobian_u(nodes_u, labels, beta, varargin) %getJacobian_f(nodes_u, labels, beta, h, n, cone_idx)
    geo = varargin{1}{1}{1};
    outer_idx = varargin{1}{1}{2};
    stat_idxs = varargin{1}{1}{3};
    rib_d = varargin{1}{1}{4};
    rib_n = varargin{1}{1}{5};
    num_outer_nodes = length(outer_idx);
    c = geo.c;
    n = geo.n;
    
    second_major_valley_indexes = find(labels < 0)';
    num_rotational_constrs = numel(second_major_valley_indexes);
    
        J = sparse(num_outer_nodes + (length(nodes_u)/3 - num_rotational_constrs) -n*(rib_n-1)*(rib_d>0)...
        +num_rotational_constrs*3 + n*(rib_d>0) + 2, length(nodes_u)); % the +2 is from fixed inner node
    
    % g_outer
    x_coord_indexes = outer_idx*3-2;
    y_coord_indexes = outer_idx*3-1;
    J(sub2ind(size(J),1:num_outer_nodes, x_coord_indexes)) = 2*nodes_u(x_coord_indexes); %2*dfdx
    J(sub2ind(size(J),1:num_outer_nodes, y_coord_indexes)) = 2*nodes_u(y_coord_indexes); %2*dfdy
    
    
    % g_rot
    major_valley_indexes = arrayfun(@(x) find(labels == ...
        -labels(x, 1), 1), second_major_valley_indexes);
    
    rotational_constr_indexes = num_outer_nodes + (1:3*num_rotational_constrs);
    Jrot = getRotationalConstraintJacobian(major_valley_indexes, ...
        second_major_valley_indexes, nodes_u, beta);
    J(rotational_constr_indexes, :) = Jrot;
    
    
    % g_surface - Apply to all but second major valley nodes
    surface_node_indexes = find(~(labels<=0))';
    if rib_d > 0
        surface_node_indexes((end-(rib_n-1)*n+1):(end-n)) = [];
    end
    surface_constraint_indexes = num_rotational_constrs * 3 + num_outer_nodes + (1:length(surface_node_indexes));
    J(surface_constraint_indexes, :) = getSurfaceConstraintJacobian(surface_node_indexes, nodes_u, c);
    
    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    J(end-1, x_coord_indexes) = 1;
    J(end, y_coord_indexes) = 1;
end