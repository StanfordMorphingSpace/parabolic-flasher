function J = getJacobian_f(nodes_u, labels, beta, varargin)
% For Jacobian F, the inputs are nodes u, correct?
h = varargin{1}{1}{1};
n = varargin{1}{1}{2};
% Do we need this variable?
cone_idx = varargin{1}{1}{3};
c = varargin{1}{1}{5};
inner_fixed_node_indices = varargin{1}{1}{6};
rib_n = varargin{1}{1}{7};
if rib_n>0
    inner_fixed_node_indices((end-rib_n+2):end) = [];
end
% Also why not constr valley nodes?
inner_fixed_node_indices = inner_fixed_node_indices(isnan(labels(inner_fixed_node_indices)));

second_major_valley_indexes = find(labels < 0)';
major_valley_indexes = arrayfun(@(x) find(labels == ...
    -labels(x, 1), 1), second_major_valley_indexes);

num_surface_constrs = length(inner_fixed_node_indices);
num_spiral_constraints = length(major_valley_indexes);
num_rotational_constrs = numel(second_major_valley_indexes);
% g_spiral
dfdx = @(x, y) y./(x.^2+y.^2) + beta*x./(2*h*sqrt(x.^2 + y.^2));
dfdy = @(x, y) -x./(x.^2+y.^2) + beta*y./(2*h*sqrt(x.^2 + y.^2));

% num_constrs = num(g_spiral)+num(g_inner)+num(g_rot)
J = sparse(num_rotational_constrs+num_spiral_constraints+num_surface_constrs, length(nodes_u));

% why not constrain all major fold line indices?
spiral_constr_indexes = 1:num_spiral_constraints;
x_coord_indexes = major_valley_indexes*3-2;
y_coord_indexes = major_valley_indexes*3-1;
J(sub2ind(size(J), spiral_constr_indexes, x_coord_indexes)) = dfdx(nodes_u(x_coord_indexes), nodes_u(y_coord_indexes));
J(sub2ind(size(J), spiral_constr_indexes, y_coord_indexes)) = dfdy(nodes_u(x_coord_indexes), nodes_u(y_coord_indexes));

% g_rot
rotational_constr_indexes = num_spiral_constraints + (1:3*num_rotational_constrs);
Jrot = getRotationalConstraintJacobian(major_valley_indexes, ...
    second_major_valley_indexes, nodes_u, beta);
J(rotational_constr_indexes, :) = Jrot;

% g_inner (applied as a surface constraint)
% Why? Also why not consider second major valley constrs?
surface_constraint_indexes = num_spiral_constraints + 3*num_rotational_constrs  + (1:num_surface_constrs);
J(surface_constraint_indexes, :) = getSurfaceConstraintJacobian(inner_fixed_node_indices, nodes_u, c);
end