function J = getJacobian_u_fast(nodes_u, labels, beta, varargin) %getJacobian_f(nodes_u, labels, beta, h, n, cone_idx)
outer_idx = varargin{1}{1}{1};
num_outer_nodes = length(outer_idx);
c = varargin{1}{1}{2};


second_major_valley_indexes = find(labels < 0)';
num_rotational_constrs = numel(second_major_valley_indexes);


% num_constrs = num(g_outer)+num(g_rot)+num(g_surface)
J = sparse(num_outer_nodes + (length(nodes_u)/3 - num_rotational_constrs) ...
    +num_rotational_constrs*3, length(nodes_u));

% g_outer
x_coord_indexes = outer_idx*3-2;
y_coord_indexes = outer_idx*3-1;
J(sub2ind(size(J),1:num_outer_nodes, x_coord_indexes)) = 2*nodes_u(x_coord_indexes); %2*dfdx
J(sub2ind(size(J),1:num_outer_nodes, y_coord_indexes)) = 2*nodes_u(y_coord_indexes); %2*dfdy


% g_rot
major_valley_indexes = arrayfun(@(x) find(labels == ...
    -labels(x, 1), 1), second_major_valley_indexes);

x_constr_indexes = num_outer_nodes + 1 +(0:3:3*num_rotational_constrs-3);
y_constr_indexes = x_constr_indexes+1;
z_constr_indexes = x_constr_indexes+2;

J(sub2ind(size(J), x_constr_indexes, major_valley_indexes*3 - 2)) = -cos(beta);
J(sub2ind(size(J), x_constr_indexes, major_valley_indexes*3 - 1)) =  sin(beta);
J(sub2ind(size(J), x_constr_indexes, 3*second_major_valley_indexes-2)) =  1;

J(sub2ind(size(J), y_constr_indexes, major_valley_indexes*3 - 2)) = -sin(beta);
J(sub2ind(size(J), y_constr_indexes, major_valley_indexes*3 - 1)) = -cos(beta);
J(sub2ind(size(J), y_constr_indexes, 3*second_major_valley_indexes-1))  =  1;

J(sub2ind(size(J), z_constr_indexes, major_valley_indexes*3))     = -1;
J(sub2ind(size(J), z_constr_indexes, 3*second_major_valley_indexes))  =  1;

% g_surface - Apply to all but second major valley nodes
node_indexes = find(~(labels<0))';
x_coord_indexes = node_indexes*3-2;
y_coord_indexes = node_indexes*3-1;
z_coord_indexes = node_indexes*3;
surface_constraint_indexes = num_rotational_constrs * 3 + num_outer_nodes + (1:length(node_indexes));
J(sub2ind(size(J), surface_constraint_indexes, x_coord_indexes)) = (2*c)*nodes_u(x_coord_indexes);
J(sub2ind(size(J), surface_constraint_indexes, y_coord_indexes)) = (2*c)*nodes_u(y_coord_indexes);
J(sub2ind(size(J), surface_constraint_indexes, z_coord_indexes)) = -1;
end