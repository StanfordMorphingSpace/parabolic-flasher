function [p_u, p_f, v_u, v_f, a_u, a_f, E_crease, E_axial, E_v_u, E_v_f] = DR_Step_parabolic(p_u, p_f, v_u, v_f, ref, labels, edges, tria1, dt, geo, mass, i, stat_idxs, outer_idx, outer_R, cone_idx)
    beta    = 2*pi()/geo.N;

    angles_u = foldedCreaseAngles(ref', p_u', edges, tria1);
    lengths_u     = getEdgeLengths(p_u(:, 1:3)', edges);

    lengths_f   = getEdgeLengths(p_f', edges);

    % unfolded
    [F_axial_u, ~, ~, ~, F_damping_u] = getDRForces(p_u, ref, edges, tria1, geo.k_axial, geo.k_fold, lengths_f, zeros(length(lengths_f)), geo.gamma, v_u, mass);

    a_u = (F_axial_u + F_damping_u)./mass;

    % folded
    [F_axial_f, E_axial, F_crease_f, E_crease, F_damping_f] = getDRForces(p_f, ref, edges, tria1, geo.k_axial, geo.k_fold, lengths_u, angles_u, geo.gamma, v_f, mass);

    a_f = (F_axial_f + F_damping_f + F_crease_f)./mass;

    for j = 1:length(stat_idxs) % exclude inner points so they stay fixed
        if j > 2
            stat_a = (a_f(stat_idxs(j), :) + a_u(stat_idxs(j), :))./2;
        else
            stat_a = [0,0,0];
        end
        a_u(stat_idxs(j), :) = stat_a; %stat_a; 
        a_f(stat_idxs(j), :) = stat_a; %stat_a;
    end

    [p_u, v_u, E_v_u] = makeStep(p_u, v_u, a_u, @getJacobian_u_fast, @getb_u, labels, beta, i, mass, dt, {outer_idx, geo.c, outer_R});
    [p_f, v_f, E_v_f] = makeStep(p_f, v_f, a_f, @getJacobian_f_fast, @getb_f, labels, beta, i, mass, dt, {geo.h, geo.n, cone_idx, geo.A, geo.c, stat_idxs});
    
end

function b = getb_u(nodes_f, labels, beta, varargin)
    outer_idx = varargin{1}{1}{1};
    c = varargin{1}{1}{2};
    outer_R = varargin{1}{1}{3};
    % constraint residuals
    % here nodes_f is a stacked column vector of all node coordinates

    b = zeros(1, 1);
    j = 1;

    for i = 1:length(outer_idx)
        idx1 = outer_idx(i)*3-2;
        idx2 = outer_idx(i)*3-1;
        
        b(j, 1) = - outer_R^2 + (nodes_f(idx1)^2 + nodes_f(idx2)^2);
        j = j+1;
    end
    

    for i = 1:3:length(nodes_f)
        % rotational symmetry
        if labels((i+2)/3, 1) < 0 % rot valley
            orig_idx = find(labels==-labels((i+2)/3, 1));
            b(j, 1) = nodes_f(i)-nodes_f(orig_idx*3-2)*cos(beta) + nodes_f(orig_idx*3-1)*sin(beta);

            b(j+1, 1) = nodes_f(i+1)-nodes_f(orig_idx*3-2)*sin(beta)-nodes_f(orig_idx*3-1)*cos(beta);

            b(j+2, 1) = -nodes_f(orig_idx*3)+nodes_f(i+2);
            j = j+3;
        else % surface constraints
            
            b(j, 1) = c*(nodes_f(i)^2 + nodes_f(i+1)^2) - nodes_f(i+2);
            j = j+1;
            
        end
    end

end

function b = getb_f(nodes_f, labels, beta, varargin)
    h = varargin{1}{1}{1};
    n = varargin{1}{1}{2};
    cone_idx = varargin{1}{1}{3};
    A = varargin{1}{1}{4};
    c = varargin{1}{1}{5};
    stat_idxs = varargin{1}{1}{6};
    % constraint residuals
    % here nodes_f is a stacked column vector of all node coordinates
    fv = @(x, y) mod(beta*(sqrt(x^2+y^2)-A)/(2*h), 2*pi) - mod(atan2(y, x), 2*pi); % valley
    fm = @(x, y) mod(beta*(sqrt(x^2+y^2)-A+h)/(2*h), 2*pi) - mod(atan2(y, x), 2*pi); % mountain

    b = zeros(1, 1);
    j = 1;
    
    % spiral
    for i = 1:3:length(nodes_f)
        if ~isnan(labels((i+2)/3, 1)) % major fold line
            if labels((i+2)/3, 1) <= n && labels((i+2)/3, 1) > 0 % valley
                b(j, 1) = fv(nodes_f(i), nodes_f(i+1));
                j = j + 1;
            % elseif labels((i+2)/3, 1) > n + cone_idx % mountain
            %     b(j, 1) = fm(nodes_f(i), nodes_f(i+1));
            %     j = j + 1;
            end

        end
    end
    % rotational symmetry
    for i = 1:3:length(nodes_f)
        if labels((i+2)/3, 1) < 0 % rot valley
            orig_idx = find(labels==-labels((i+2)/3, 1));
            b(j, 1) = nodes_f(i)-nodes_f(orig_idx*3-2)*cos(beta) + nodes_f(orig_idx*3-1)*sin(beta);

            b(j+1, 1) = nodes_f(i+1)-nodes_f(orig_idx*3-2)*sin(beta)-nodes_f(orig_idx*3-1)*cos(beta);

            b(j+2, 1) = -nodes_f(orig_idx*3)+nodes_f(i+2);
            j = j+3;
        else % center surface constraints
            if ismember((i+2)/3, stat_idxs(3:end))
                b(j, 1) = c*(nodes_f(i)^2 + nodes_f(i+1)^2) - nodes_f(i+2);
                j = j+1;
            end
        end
    end
end