function b = getSurfaceConstraints(node_indexes, node_coords, surf_func, offset)
x_coord_indexes = node_indexes*3-2;
y_coord_indexes = x_coord_indexes+1;
z_coord_indexes = x_coord_indexes+2;
r = sqrt(node_coords(x_coord_indexes).^2 + node_coords(y_coord_indexes).^2);
b = surf_func(r) + offset - node_coords(z_coord_indexes);
end