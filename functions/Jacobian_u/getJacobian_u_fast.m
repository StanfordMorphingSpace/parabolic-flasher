function J = getJacobian_u_fast(nodes_u, labels, beta, varargin) %getJacobian_f(nodes_u, labels, beta, h, n, cone_idx)
    outer_idx = varargin{1}{1}{1};
    num_outer_nodes = length(outer_idx);
    c = varargin{1}{1}{2};

    % jacobian of equality constraints (g_outer)
    % num_constraints = num(g_outer)+num(g_surface)+num(g_rot)
    J = sparse(num_outer_nodes + length(nodes_u)/3, length(nodes_u));
    x_coord_indexes = outer_idx*3-2;
    y_coord_indexes = outer_idx*3-1;
    J(sub2ind(size(J),1:num_outer_nodes, x_coord_indexes)) = 2*nodes_u(x_coord_indexes); %2*dfdx
    J(sub2ind(size(J),1:num_outer_nodes, y_coord_indexes)) = 2*nodes_u(y_coord_indexes); %2*dfdy
    
    j = num_outer_nodes+1;
    for i = 1:3:length(nodes_u)
        % rotational symmetry
        if labels((i+2)/3, 1) < 0 % rot valley
            orig_idx = find(labels==-labels((i+2)/3, 1));
            J(j, orig_idx*3-2) = -cos(beta);
            J(j, orig_idx*3-1) = sin(beta);
            J(j, i) = 1;

            J(j+1, orig_idx*3-2) = -sin(beta);
            J(j+1, orig_idx*3-1) = -cos(beta);
            J(j+1, i+1) = 1;

            J(j+2, orig_idx*3) = -1;
            J(j+2, i+2) = 1;
            j = j+3;
        else 
            % surface constraints
            J(j, i) = 2*c*nodes_u(i);
            J(j, i+1) = 2*c*nodes_u(i+1);
            J(j, i+2) = -1;
            j = j+1;
        end
    end
   
end