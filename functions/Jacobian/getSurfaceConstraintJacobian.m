function J = getSurfaceConstraintJacobian(node_indexes, node_coords, surf_func_prime)
num_constrs = length(node_indexes);
num_vars = length(node_coords);
x_coord_indexes = node_indexes*3-2;
y_coord_indexes = x_coord_indexes+1;
z_coord_indexes = x_coord_indexes+2;
constraint_indexes = 1:num_constrs;

x = node_coords(x_coord_indexes)';
y = node_coords(y_coord_indexes)';
r = sqrt(x.^2 + y.^2);
dzdr = surf_func_prime(r);

J = sparse([constraint_indexes, constraint_indexes, constraint_indexes], ...
    [x_coord_indexes, y_coord_indexes, z_coord_indexes], ...
    [dzdr.*x./r, dzdr.*y./r, -1*ones(1, num_constrs)], ...
    num_constrs, num_vars);
end