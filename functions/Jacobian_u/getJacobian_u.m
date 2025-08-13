function J = getJacobian_u(nodes_f, labels, beta, varargin) %getJacobian_f(nodes_f, labels, beta, h, n, cone_idx)
    outer_idx = varargin{1}{1}{1};
    c = varargin{1}{1}{2};
    % jacobian of equality constraints
    % here nodes_f is a stacked column vector of all node coordinates
    dfdx = @(x, y) 2*x;
    dfdy = @(x, y) 2*y;
    num_free_nodes = length(outer_idx);

    J = zeros(num_free_nodes, length(nodes_f));
    j = 1;

    for i = 1:length(outer_idx)
        idx1 = outer_idx(i)*3-2;
        idx2 = outer_idx(i)*3-1;
        
        J(j, idx1) = dfdx(nodes_f(idx1), nodes_f(idx2));
        J(j, idx2) = dfdy(nodes_f(idx1), nodes_f(idx2));
        j = j+1;
    end
    
    for i = 1:3:length(nodes_f)
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
        else % surface constraints
            J(j, i) = 2*c*nodes_f(i);
            J(j, i+1) = 2*c*nodes_f(i+1);
            J(j, i+2) = -1;
            j = j+1;
        end
    end
   
end