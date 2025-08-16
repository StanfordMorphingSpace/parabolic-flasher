function J = getJacobian_u(nodes_f, labels, beta, varargin) %getJacobian_f(nodes_f, labels, beta, h, n, cone_idx)
outer_idx = varargin{1}{1}{1};
num_outer_nodes = length(outer_idx);
c = varargin{1}{1}{2};


second_major_valley_indexes = find(labels < 0)';
num_rotational_constrs = numel(second_major_valley_indexes);


% num_constrs = num(g_outer)+num(g_rot)+num(g_surface)
J = sparse(num_outer_nodes + (length(nodes_f)/3 - num_rotational_constrs) ...
    +num_rotational_constrs*3, length(nodes_f));

% g_outer
x_coord_indexes = outer_idx*3-2;
y_coord_indexes = outer_idx*3-1;
J(sub2ind(size(J),1:num_outer_nodes, x_coord_indexes)) = 2*nodes_f(x_coord_indexes); %2*dfdx
J(sub2ind(size(J),1:num_outer_nodes, y_coord_indexes)) = 2*nodes_f(y_coord_indexes); %2*dfdy


% g_rot
major_valley_indexes = arrayfun(@(x) find(labels == ...
    -labels(x, 1), 1), second_major_valley_indexes);

rotational_constr_indexes = num_outer_nodes + (1:3*num_rotational_constrs);
Jrot = getRotationalConstraintJacobian(major_valley_indexes, ...
    second_major_valley_indexes, nodes_f, beta);
J(rotational_constr_indexes, :) = Jrot;


% g_surface - Apply to all but second major valley nodes
surface_node_indexes = find(~(labels<0))';
surface_constraint_indexes = num_rotational_constrs * 3 + num_outer_nodes + (1:length(surface_node_indexes));
J(surface_constraint_indexes, :) = getSurfaceConstraintJacobian(surface_node_indexes, nodes_f, c);
end