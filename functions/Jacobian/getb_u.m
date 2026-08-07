
function b = getb_u(nodes_u, labels, beta, varargin)
    geo = varargin{1}{1}{1};
    outer_idx = varargin{1}{1}{2};
    stat_idxs = varargin{1}{1}{3};
    rib_d = varargin{1}{1}{4};
    rib_n = varargin{1}{1}{5};
    num_outer_nodes = length(outer_idx);
    c = geo.c;
    n = geo.n;
    A = geo.A;
    outer_R = geo.R;
    
    second_major_valley_indexes = find(labels < 0)';
    num_rotational_constrs = numel(second_major_valley_indexes);
    
    b = zeros(num_outer_nodes + (length(nodes_u)/3 - num_rotational_constrs) -n*(rib_n-1)*(rib_d>0)...
        +num_rotational_constrs*3 + n*(rib_d>0) +2, 1);
    
    % g_outer
    x_coord_indexes = outer_idx.*3-2;
    y_coord_indexes = x_coord_indexes+1;
    b(1:num_outer_nodes)= -outer_R^2 + (nodes_u(x_coord_indexes).^2 + nodes_u(y_coord_indexes).^2);
    
    % g_rot
    major_valley_indexes = arrayfun(@(x) find(labels == ...
        -labels(x, 1), 1), second_major_valley_indexes);
    num_rotational_constrs = length(second_major_valley_indexes);
    b(num_outer_nodes+(1:3*num_rotational_constrs)) = getRotationalConstraints(major_valley_indexes, second_major_valley_indexes, nodes_u, beta);
    
    %g_surface
    surface_node_indexes = find(~(labels<=0))';
    if rib_d >0
        surface_node_indexes((end-(rib_n-1)*n+1):end) = [];
    end
    surface_constraint_indexes = num_rotational_constrs * 3 + num_outer_nodes + (1:length(surface_node_indexes));
    b(surface_constraint_indexes)=  getSurfaceConstraints(surface_node_indexes, nodes_u, c, 0);
    
    if rib_d > 0
        rib_edge_indexes = (length(nodes_u)/3) - 0:(n-1);
        rib_constraint_indexes = num_rotational_constrs * 3 + num_outer_nodes + length(surface_node_indexes) + 1:n;
        b(rib_constraint_indexes)=  getSurfaceConstraints(rib_edge_indexes, nodes_u, c, -rib_d);
    end

    % constrain inner node
    x_coord_indexes = stat_idxs(1)*3-2;
    y_coord_indexes = x_coord_indexes+1;
    b(end-1) = nodes_u(x_coord_indexes) - A;
    b(end) = nodes_u(y_coord_indexes);

end
