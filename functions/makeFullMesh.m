function [unfolded, edges_one, faces_one] = makeFullMesh(nodes_u, edges, faces, rot, N)

    unfolded = nodes_u;
    faces_one = faces;
    
    for j = 1:(N-1)
        rot_u = (rot^j*nodes_u')';
    
        tria_temp = zeros(size(faces));
        
        for i = 1:length(rot_u)
            found_idx = find(ismembertol(unfolded, rot_u(i, :), 1e-4, 'ByRows', true));
            tria_idx = find(faces(:, 1:3)==i);
            if isempty(found_idx)
                unfolded = [unfolded; rot_u(i, :)];
                tria_temp(tria_idx) = length(unfolded);
            else
                tria_temp(tria_idx) = found_idx;
            end
            
        end
        
        faces_one = [faces_one; tria_temp];
    end
        
    % assign edges
    for i = 1:length(faces_one)
        edge1 = faces_one(i, 1:2);
        edge2 = faces_one(i, 2:3);
        edge3 = faces_one(i, [1 3]);
        if i == 1
            edges_one = [edge1; edge2; edge3];
        else
            if isempty(find(ismember(edges_one, edge1, 'rows'))) & isempty(find(ismember(edges_one, flip(edge1), 'rows')))
                edges_one = [edges_one; edge1];
            end
            if isempty(find(ismember(edges_one, edge2, 'rows'))) & isempty(find(ismember(edges_one, flip(edge2), 'rows')))
                edges_one = [edges_one; edge2];
            end
            if isempty(find(ismember(edges_one, edge3, 'rows'))) & isempty(find(ismember(edges_one, flip(edge3), 'rows')))
                edges_one = [edges_one; edge3];
            end
        end
    end

end


