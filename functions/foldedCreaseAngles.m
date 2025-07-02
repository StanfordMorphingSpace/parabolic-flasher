function angles = foldedCreaseAngles_correct(nodes_deployed, nodes_folded, edges, faces)
    nEdges = size(edges, 1);    
    angles = nan(1, nEdges);
    
    for i = 1:nEdges
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);
        
        faces_with_node1        = (faces(:, 1) == p1_index) | (faces(:, 2) == p1_index) | (faces(:, 3) == p1_index);
        faces_with_node2        = (faces(:, 1) == p2_index) | (faces(:, 2) == p2_index) | (faces(:, 3) == p2_index);
        faces_adjacent_to_edge  = faces((faces_with_node1 & faces_with_node2), :);
        
        if size(faces_adjacent_to_edge, 1) ~= 2
            % disp([num2str(size(faces_adjacent_to_edge, 1)) ' faces at edge connecting nodes ' num2str(p1_index) ' and ' num2str(p2_index)])            
            angles(i) = nan;
        else
            face1 = faces_adjacent_to_edge(1, :);
            face2 = faces_adjacent_to_edge(2, :);
            
            p3_1_index = face1((face1 ~= p1_index) & (face1 ~= p2_index));
            p3_2_index = face2((face2 ~= p1_index) & (face2 ~= p2_index));

            p1 = nodes_folded(1:3, p1_index);
            p2 = nodes_folded(1:3, p2_index);
            p1_u = nodes_deployed(1:3, p1_index);
            p2_u = nodes_deployed(1:3, p2_index);
                        
            p3_1 = nodes_folded(1:3, p3_1_index);
            p3_2 = nodes_folded(1:3, p3_2_index); 
            p3_1_u = nodes_deployed(1:3, p3_1_index);
            p3_2_u = nodes_deployed(1:3, p3_2_index);

            n1_u = cross(p1_u-p3_1_u, p2_u-p3_1_u); n1_u=n1_u/norm(n1_u);
            n2_u = cross(p1_u-p3_2_u, p2_u-p3_2_u); n2_u=n2_u/norm(n2_u);

            n1 = cross(p1-p3_1, p2-p3_1); n1=n1/norm(n1);
            n2 = cross(p3_2-p1, p3_2-p2); n2=n2/norm(n2);

            if n1_u(3) < 0
                n1 = -n1;
            end

            if n2_u(3) < 0
                n2 = -n2;
            end

            angles(i) = angleBetweenVectors3d(n1, n2); angles(i) = real(angles(i));

            if dot(n1, p3_2-p3_1) < 0
                angles(i) = -angles(i);
            end

            if abs(angles(i)) == pi
                angles(i) = 0;
            end
        end
    end
end