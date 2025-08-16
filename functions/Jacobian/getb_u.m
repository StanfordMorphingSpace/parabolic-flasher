
function b = getb_u(nodes_u, labels, beta, varargin)
outer_idx = varargin{1}{1}{1};
c = varargin{1}{1}{2};
outer_R = varargin{1}{1}{3};
num_outer_nodes = length(outer_idx);

second_major_valley_indexes = find(labels < 0)';
num_rotational_constrs = numel(second_major_valley_indexes);

% num_constrs = num(g_outer)+num(g_rot)+num(g_surface)
b = zeros(num_outer_nodes + (length(nodes_u)/3 - num_rotational_constrs) ...
    + num_rotational_constrs*3, 1);


% g_outer
x_coord_indexes = outer_idx.*3-2;
y_coord_indexes = x_coord_indexes+1;
b(1:num_outer_nodes)= -outer_R^2 + (nodes_u(x_coord_indexes).^2 + nodes_u(y_coord_indexes).^2);

% g_rot
major_valley_indexes = arrayfun(@(x) find(labels == ...
    -labels(x, 1), 1), second_major_valley_indexes);
num_rotational_constrs = length(second_major_valley_indexes);
b(1:3*num_rotational_constrs) = getRotationalConstraints(major_valley_indexes, second_major_valley_indexes, nodes_u, beta);

%g_surface
surface_node_indexes = find(~(labels<0))';
surface_constraint_indexes = num_rotational_constrs * 3 + num_outer_nodes + (1:length(surface_node_indexes));
b(surface_constraint_indexes)=  getSurfaceConstraints(surface_node_indexes, nodes_u, c);
end
