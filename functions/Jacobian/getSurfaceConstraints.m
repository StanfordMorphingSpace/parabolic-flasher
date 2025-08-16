function b = getSurfaceConstraints(node_indexes, node_coords, c)
x_coord_indexes = node_indexes*3-2;
y_coord_indexes = x_coord_indexes+1;
z_coord_indexes = x_coord_indexes+2;
b = c*(node_coords(x_coord_indexes).^2 + node_coords(y_coord_indexes).^2) - node_coords(z_coord_indexes);
end