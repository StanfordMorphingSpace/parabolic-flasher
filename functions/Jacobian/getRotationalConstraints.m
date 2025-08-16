function b = getRotationalConstraints(major_valley_indexes, second_major_valley_indexes, node_coords, beta)
num_rotational_constrs = length(second_major_valley_indexes);
b=zeros(num_rotational_constrs*3,1);
x_constr_indexes = (1:3:3*num_rotational_constrs);
y_constr_indexes = x_constr_indexes+1;
z_constr_indexes = x_constr_indexes+2;
second_major_valley_x_coord_idxs = second_major_valley_indexes*3-2;
second_major_valley_y_coord_idxs = second_major_valley_x_coord_idxs+1;
second_major_valley_z_coord_idxs = second_major_valley_y_coord_idxs+1;

major_valley_x_coord_idxs = major_valley_indexes*3-2;
major_valley_y_coord_idxs = major_valley_x_coord_idxs+1;
major_valley_z_coord_idxs = major_valley_y_coord_idxs+1;

b(x_constr_indexes) = node_coords(second_major_valley_x_coord_idxs) - node_coords(major_valley_x_coord_idxs)*cos(beta) + node_coords(major_valley_y_coord_idxs)*sin(beta);
b(y_constr_indexes) = node_coords(second_major_valley_y_coord_idxs) - node_coords(major_valley_x_coord_idxs)*sin(beta) - node_coords(major_valley_y_coord_idxs)*cos(beta);
b(z_constr_indexes) = node_coords(second_major_valley_z_coord_idxs) - node_coords(major_valley_z_coord_idxs);
end