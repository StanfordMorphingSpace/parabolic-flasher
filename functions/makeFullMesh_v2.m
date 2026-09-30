function [nodes_one, nodes_f_one, edges_one, faces_one, angles_one] = makeFullMesh_v2(nodes_u, nodes_f, edges, faces, angles, rot, N)
    
    % Ensure angles is a column vector
    angles = angles(:);
    
    % Number of nodes, edges, and faces in a single sector
    num_nodes = size(nodes_u, 1);
    num_edges = size(edges, 1);
    num_faces = size(faces, 1);
    
    % Preallocate global arrays for speed
    all_nodes_u = zeros(num_nodes * N, 3);
    all_nodes_f = zeros(num_nodes * N, 3);
    all_faces = zeros(num_faces * N, 3);
    all_edges = zeros(num_edges * N, 2);
    all_angles = repmat(angles, N, 1);
    
    % 1. Generate all rotated sectors for both unfolded and folded states
    current_rot = eye(3);
    for j = 1:N
        node_idx = (j-1)*num_nodes + 1 : j*num_nodes;
        face_idx = (j-1)*num_faces + 1 : j*num_faces;
        edge_idx = (j-1)*num_edges + 1 : j*num_edges;
        
        % Rotate both unfolded and folded nodes identically
        all_nodes_u(node_idx, :) = (current_rot * nodes_u')';
        all_nodes_f(node_idx, :) = (current_rot * nodes_f')';
        
        % Offset faces and edges to match the current sector block
        all_faces(face_idx, :) = faces + (j-1)*num_nodes;
        all_edges(edge_idx, :) = edges + (j-1)*num_nodes;
        
        current_rot = current_rot * rot;
    end
    
    % 2. Find unique nodes while FORCING stable original order
    % uniquetol sorts by default, so we extract the indices to reconstruct the stable order
    [~, ia, ic_sorted] = uniquetol(all_nodes_u, 1e-4, 'ByRows', true, 'DataScale', 1);
    
    % Sort 'ia' (index of first occurrences) to restore sequential appending order
    [ia_stable, stable_idx] = sort(ia);
    
    % Extract unique nodes in their original, un-jumbled order
    nodes_one = all_nodes_u(ia_stable, :);
    nodes_f_one = all_nodes_f(ia_stable, :); % Exact same slicing applied to folded nodes
    
    % Reconstruct the inverse mapping (ic) to match the new stable order
    remap = zeros(length(ia), 1);
    remap(stable_idx) = 1:length(ia);
    ic_nodes = remap(ic_sorted);
    
    % 3. Remap face and edge connectivity cleanly
    faces_one = ic_nodes(all_faces);
    mapped_edges = ic_nodes(all_edges);
    
    % 4. Standardize edge directions and extract unique edges (STABLE)
    mapped_edges = sort(mapped_edges, 2);
    [edges_one, ia_edges, ic_edges] = unique(mapped_edges, 'rows', 'stable');
    angles_one = all_angles(ia_edges);
    
    % 5. Handle overlapping edge angles
    % Count occurrences of each edge in the stable index mapping
    edge_counts = accumarray(ic_edges, 1);
    angles_one(edge_counts > 1) = 0;
    
    % 6. Clean up duplicate faces (STABLE)
    % Sort rows temporarily just to find duplicates, but preserve original face order
    sorted_faces = sort(faces_one, 2);
    [~, unique_face_idx, ~] = unique(sorted_faces, 'rows', 'stable');
    faces_one = faces_one(unique_face_idx, :);

end