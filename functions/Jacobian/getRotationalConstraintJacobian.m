function J = getRotationalConstraintJacobian(major_valley_indexes, second_major_valley_indexes, node_coords, beta)
num_rotational_constrs = numel(second_major_valley_indexes);
J = sparse(3*num_rotational_constrs, length(node_coords));

x_constr_indexes = 1:3:3*num_rotational_constrs;
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
end