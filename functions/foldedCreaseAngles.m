function angles = foldedCreaseAngles(nodes_deployed, nodes_folded, edges, adj_faces)
    nEdges = size(edges, 1);    
    angles = nan(nEdges, 1);
    
    for i = 1:nEdges
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);
        
        if length(adj_faces{i}) ~= 2
            angles(i) = nan;
        else            
            p3_1_index = adj_faces{i}(1);
            p3_2_index = adj_faces{i}(2);

            p1 = nodes_folded(p1_index, :);
            p2 = nodes_folded(p2_index, :);
            p1_u = nodes_deployed(p1_index, :);
            p2_u = nodes_deployed(p2_index, :);
                        
            p3_1 = nodes_folded(p3_1_index, :);
            p3_2 = nodes_folded(p3_2_index, :); 
            p3_1_u = nodes_deployed(p3_1_index, :);
            p3_2_u = nodes_deployed(p3_2_index, :);

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