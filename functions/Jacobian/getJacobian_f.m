function J = getJacobian_f(nodes_f, labels, beta, varargin)
h = varargin{1}{1}{1};
n = varargin{1}{1}{2};
cone_idx = varargin{1}{1}{3};
c = varargin{1}{1}{5};
stat_idxs = varargin{1}{1}{6};
% jacobian of equality constraints
% here nodes_f is a stacked column vector of all node coordinates
dfdx = @(x, y) y/(x^2+y^2) + beta*x/(2*h*sqrt(x^2 + y^2));
dfdy = @(x, y) -x/(x^2+y^2) + beta*y/(2*h*sqrt(x^2 + y^2));

J = zeros(1, length(nodes_f));
j = 1;

% spiral constraints
for i = 1:3:length(nodes_f)
    if ~isnan(labels((i+2)/3, 1)) % major fold line
        if (labels((i+2)/3, 1) <= n && labels((i+2)/3, 1) > 0) % || labels((i+2)/3, 1) > n + cone_idx % everything else
            J(j, i) = dfdx(nodes_f(i), nodes_f(i+1));
            J(j, i+1) = dfdy(nodes_f(i), nodes_f(i+1));
            j = j + 1;
        end
    end
end
% rotational symmetry
for i = 1:3:length(nodes_f)
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
    end
end
for i = 1:3:length(nodes_f)
    if ~(labels((i+2)/3, 1) < 0) % rot valley
        % surface constraints
        if ismember((i+2)/3, stat_idxs(3:end))
            J(j, i) = 2*c*nodes_f(i);
            J(j, i+1) = 2*c*nodes_f(i+1);
            J(j, i+2) = -1;
            j = j+1;
        end
    end
end
end
