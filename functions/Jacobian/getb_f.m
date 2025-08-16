function b = getb_f(nodes_f, labels, beta, varargin)
h = varargin{1}{1}{1};
n = varargin{1}{1}{2};
cone_idx = varargin{1}{1}{3};
A = varargin{1}{1}{4};
c = varargin{1}{1}{5};
inner_fixed_node_indices = varargin{1}{1}{6};
inner_fixed_node_indices = inner_fixed_node_indices(isnan(labels(inner_fixed_node_indices)));

second_major_valley_indexes = find(labels < 0)';
major_valley_indexes = arrayfun(@(x) find(labels == ...
    -labels(x, 1), 1), second_major_valley_indexes);
sorted_major_valley_indexes = sort(major_valley_indexes);

num_surface_constrs = length(inner_fixed_node_indices);
num_spiral_constraints = length(major_valley_indexes);
num_rotational_constrs = numel(second_major_valley_indexes);
% num_constrs = num(g_spiral)+num(g_inner)+num(g_rot)
b = ones(num_rotational_constrs+num_spiral_constraints+num_surface_constrs,1);

% g_spiral

fv = @(x, y) mod(beta*(sqrt(x.^2+y.^2)-A)./(2*h), 2*pi) - mod(atan2(y, x), 2*pi); % valley

spiral_constr_indexes = 1:num_spiral_constraints;
x_coord_indexes = sorted_major_valley_indexes*3-2;
y_coord_indexes = x_coord_indexes+1;
b(spiral_constr_indexes) = fv(nodes_f(x_coord_indexes), nodes_f(y_coord_indexes));

% g_rot
b(1:3*num_rotational_constrs) = getRotationalConstraints(major_valley_indexes, second_major_valley_indexes, nodes_f, beta);

% g_inner (Applied as g_surface)
inner_fixed_node_indices = sort(inner_fixed_node_indices);
surface_constraint_indexes = num_spiral_constraints + 3*num_rotational_constrs  + (1:num_surface_constrs);
b(surface_constraint_indexes)=  getSurfaceConstraints(inner_fixed_node_indices, nodes_f, c);
end