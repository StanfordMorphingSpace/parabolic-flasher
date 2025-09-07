function [vert_u, vert_f, vert_ref, labels, edges, faces, stat_idxs, outer_idx, cone_idx] = generateMesh_ribs(A, N, h, n, R, c, plot_on, rib_d)
    % This is all a mess and a problem for later me

    beta    = 2*pi()/N;
    delta_z = 2*A*sin(beta/2)*tan(beta/2); % changes in height of mountain for a zero-thickness flasher
    
    d       = h/cos(beta/2);
    rot     = [ cos(beta), -sin(beta), 0;...
                sin(beta), cos(beta), 0;...
                0, 0, 1];
    rot_gen  = @(theta) [ cos(theta), -sin(theta), 0;...
                sin(theta), cos(theta), 0;...
                0, 0, 1];

    x_end = @(p1, p2) ((p1(2)-p2(2))*(p2(1)*p1(2) - p1(1)*p2(2)) + sqrt((p1(1)-p2(1))^2*(R^2*norm(p1-p2)^2 - (p2(1)*p1(2) - p1(1)*p2(2))^2)))/(norm(p1-p2)^2);
    y_end = @(p1, p2) (-p2(1)^3*p1(2) + p1(1)^3*p2(2) + p1(1)*p2(1)^2*(2*p1(2) + p2(2)) - p1(1)^2*p2(1)*(p1(2) + 2*p2(2)) + (p1(2)-p2(2))*sqrt((p1(1)-p2(1))^2*(R^2*norm(p1-p2)^2 - (p2(1)*p1(2) - p1(1)*p2(2))^2)))/((p1(1)-p2(1))*norm(p1-p2)^2);
    
    init_vec = [linspace(A, R, 2); zeros(1, 2); zeros(1, 2)];
    v_init = rot_gen((pi-beta)/2)*(init_vec-init_vec(:, 1))+[A;0;0];
    m_init = rot_gen(pi/2)*(init_vec-init_vec(:, 1))+[A;0;0];
    
    v_end = [x_end(v_init(:, 1), v_init(:, 2)); y_end(v_init(:, 1), v_init(:, 2)); 0];
    %m_end = [x_end(m_init(:, 1), m_init(:, 2)); y_end(m_init(:, 1), m_init(:, 2)); 0];
    m_end = [m_init(1, 1), sqrt(R^2-m_init(1, 1)^2)];
    
    v_plane = [linspace(A, v_end(1), n); linspace(0, v_end(2), n); zeros(1, n)];
    m_plane = [linspace(A, m_end(1), n); linspace(0, m_end(2), n); zeros(1, n)];
    
    % remove tip so it's not subdivided
    idx = 0;
    dis = 0;
    while dis <= (R-A)/n
        idx = idx+1;
        dis = norm(v_plane(:, idx)-m_plane(:, idx));
    end
    idx = idx-1;
    
    % gore 1
    end_edge_d = norm(v_plane(:, end)-m_plane(:, end));
    end_n = round(end_edge_d/((R-A)/n))+1;
    end_theta = angleBetweenVectors(v_plane(:, end)', m_plane(:, end)');
    end_thetas = linspace(0, end_theta, end_n);
    end_edge = [];
    for i = 2:(end_n-1)
        end_edge = [end_edge, rot_gen(-end_thetas(i))*m_plane(:, end)];
    end
    
    nodes1 = [flip(v_plane(:, idx:end), 2), m_plane(:, idx:end), end_edge];
    
    gore1 = [length(nodes1), 1];
    for i = 1:(length(nodes1)-1)
        gore1 = [gore1; i, i+1];
    end
        
    % opts.kind = 'delaunay';
    % opts.rho2 = +1 ;
    opts.ref1 = 'preserve';
    hfun = +(R-A)/n;
    
    [vert1,etri1, ...
        tria1,tnum1] = refine2(nodes1(1:2, :)',gore1, [], opts, hfun) ;
    alphas_test = real(getAlphas([vert1, zeros(length(vert1), 1)], tria1));
    bad_idxs = [];
    for i = 1:length(alphas_test)
        if abs(alphas_test(i, 1)) < 1e-6 || abs(alphas_test(i, 2)) < 1e-6 || abs(alphas_test(i, 3)) < 1e-6
            bad_idxs = [bad_idxs, i];
        end
    end
    tria1(bad_idxs, :) = [];
    
    % add back tip
    if idx > 2
        num_vert = size(vert1, 1);
        vert1 = [vert1; v_plane(1:2, (1:idx-1))'; m_plane(1:2, (2:idx-1))'];
        tip_idx1 = find(ismember(vert1, v_plane(1:2, idx)', 'rows'));
        tip_idx2 = find(ismember(vert1, m_plane(1:2, idx)', 'rows'));
        
        tria1 = [tria1; num_vert+1, num_vert+2, num_vert+idx];
        for i = 2:(idx-1)
            if i == idx-1
                tria1 = [tria1; num_vert+i, num_vert+idx+(i-2), tip_idx2;
                           num_vert+i, tip_idx2, tip_idx1];
            else
                tria1 = [tria1; num_vert+i, num_vert+idx+(i-2), num_vert+idx+(i-1);
                       num_vert+i, num_vert+idx+(i-1), num_vert+i+1];
            end
        end
    
    elseif idx == 2
        num_vert = size(vert1, 1);
        vert1 = [vert1; v_plane(1:2, 1)'];
        tip_idx1 = find(ismember(vert1, v_plane(1:2, idx)', 'rows'));
        tip_idx2 = find(ismember(vert1, m_plane(1:2, idx)', 'rows'));
    
        tria1 = [tria1; num_vert+1, tip_idx1, tip_idx2];
    end
    
    % gore 2
    end_edge_d = norm(rot*v_plane(:, 1)-m_plane(:, 1));
    div3 = linspace(0, end_edge_d, round(end_edge_d/((R-A)/n))+1);
    div3 = div3(2:(end-1));
    dir0 = m_plane(:, 1) - rot*v_plane(:, 1); dir0 = dir0/norm(dir0);
    inner_edge = rot*v_plane(:, 1) + dir0*div3;
    
    end_edge_d = norm(rot*v_plane(:, end)-m_plane(:, end));
    end_n = round(end_edge_d/((R-A)/n))+1;
    end_theta = angleBetweenVectors((rot*v_plane(:, end))', m_plane(:, end)');
    end_thetas = linspace(0, end_theta, end_n);
    end_edge = [];
    for i = 2:(end_n-1)
        end_edge = [end_edge, rot_gen(end_thetas(i))*rot*v_plane(:, end)];
    end
    
    edge_m = [flip(inner_edge, 2), rot*v_plane(:, 1:(end))];
    nodes2 = [flip(m_plane(:, idx:end), 2), edge_m(:, (idx-1):end), end_edge];
    
    gore2 = [length(nodes2), 1];
    for i = 1:(length(nodes2)-1)
        gore2 = [gore2; i, i+1];
    end
    
    %opts.kind = 'delaunay';
    %opts.rho2 = +1 ;
    opts.ref1 = 'preserve';
    hfun = +(R-A)/n;
    
    [vert2,etri2, ...
        tria2,tnum2] = refine2(nodes2(1:2, :)',gore2, [], opts, hfun) ;
    
    alphas_test = real(getAlphas([vert2, zeros(length(vert2), 1)], tria2));
    bad_idxs = [];
    for i = 1:length(alphas_test)
        if abs(alphas_test(i, 1)) < 1e-6 || abs(alphas_test(i, 2)) < 1e-6 || abs(alphas_test(i, 3)) < 1e-6
            bad_idxs = [bad_idxs, i];
        end
    end
    tria2(bad_idxs, :) = [];
    
    % add back tip
    if idx > 2
        num_vert = size(vert2, 1);
        edge_m = [flip(inner_edge, 2), rot*v_plane];
        vert2 = [vert2;  m_plane(1:2, 1:(idx-1))'; edge_m(1:2, 1:(idx-2))'];
        tip_idx2 = find(ismember(vert2, edge_m(1:2, idx-1)', 'rows'));
        tip_idx1 = find(ismember(vert2, m_plane(1:2, idx)', 'rows'));
        
        etri2 = [etri2; num_vert+1, num_vert+idx];
        for i = 1:(idx-2)
            etri2 = [etri2; num_vert+i, num_vert+i+1; num_vert+idx+i-1, num_vert+i+idx];
        end
        tria2 = [tria2; num_vert+1, num_vert+2, num_vert+idx];
        for i = 2:(idx-1)
            if i == idx-1
                tria2 = [tria2; num_vert+i, num_vert+idx+(i-2), tip_idx2;
                           num_vert+i, tip_idx2, tip_idx1];
            else
                tria2 = [tria2; num_vert+i, num_vert+idx+(i-2), num_vert+idx+(i-1);
                       num_vert+i, num_vert+idx+(i-1), num_vert+i+1];
            end
        end
    
    elseif idx == 2
        num_vert = size(vert2, 1);
        vert2 = [vert2; v_plane(1:2, 1)'];
        tip_idx2 = find(ismember(vert2, inner_edge(1:2, end-idx+2)', 'rows'));
        tip_idx1 = find(ismember(vert2, m_plane(1:2, idx)', 'rows'));
    
        etri2 = [etri2; num_vert+1, tip_idx1; num_vert+1, tip_idx2];
        tria2 = [tria2; num_vert+1, tip_idx1, tip_idx2];
    end
    
    %% assign nodes to ruling triangles
    
    % establish ruling triangles
    ruling_nodes1 = [v_plane, m_plane(:, 2:end)];
    ruling_tris1 = [1, 2, 1+n];
    for i = 1:(n-2)
        ruling_tris1 = [ruling_tris1; i+1, i+n, i+2];
        ruling_tris1 = [ruling_tris1; i+n, i+2, i+n+1];
    end
    ruling_tris1 = [ruling_tris1; i+2, i+n+1, i+3];
    
    m_plane_ext = m_plane;
    for i = 1:length(inner_edge)
        m_plane_ext = [m_plane_ext, m_plane(:, end)+i.*(m_plane(:, 2)-m_plane(:, 1))];
    end
    
    ruling_nodes2 = [m_plane_ext, flip(inner_edge, 2), rot*v_plane(:, 1:end)];
    n_ext = n+length(inner_edge);
    ruling_tris2 = [1, 1+n_ext, 2];
    i = 1;
    while max(ruling_tris2) < length(ruling_nodes2)
        ruling_tris2 = [ruling_tris2; i+n_ext, i+1, i+n_ext+1];
        ruling_tris2 = [ruling_tris2; i+1, i+n_ext+1, i+2 ];
        i = i+1;
    end
    %ruling_tris2 = [ruling_tris2; i+n_ext+1, i+2, i+n_ext+2];
    %ruling_tris2 = [ruling_tris2; i+2, i+n_ext+2, i+3 ];
    % for i = 1:(length(inner_edge)-1)
    %     ruling_tris2 = [ruling_tris2; i+2*n-1, n, i+2*n];
    % end
    
    % asign nodes to ruling triangles
    vert1 = [vert1, zeros(length(vert1), 1)];
    labels1 = zeros(length(vert1), 2);
    for i = 1:length(vert1)
        ruling_idx = find(ismember(ruling_nodes1(1:2, :)', vert1(i, 1:2), 'rows'));
        if ~isempty(ruling_idx)
            labels1(i, 2) = ruling_idx;
        else
            labels1(i, 2) = nan;
        end
        for j = 1:length(ruling_tris1)
            [in, on] = inpolygon(vert1(i, 1), vert1(i, 2), ruling_nodes1(1, ruling_tris1(j, :)), ruling_nodes1(2, ruling_tris1(j, :)));
            if in == 1
                if on == 0
                    labels1(i, 1) = j;
                elseif isempty(ruling_idx)
                    labels1(i, 1) = j;
                else
                    labels1(i, 1) = nan;
                end
                break
            end
        end
    end
    
    labels1(labels1==0) = max(labels1(:, 1)); % bug workaround
    
    vert2 = [vert2, zeros(length(vert2), 1)];
    labels2 = zeros(length(vert2), 2);
    for i = 1:length(vert2)
        ruling_idx = find(ismember(ruling_nodes2(1:2, :)', vert2(i, 1:2), 'rows'));
        if ~isempty(ruling_idx)
            labels2(i, 2) = ruling_idx;
        else
            labels2(i, 2) = nan;
        end
        for j = 1:length(ruling_tris2)
            [in, on] = inpolygon(vert2(i, 1), vert2(i, 2), ruling_nodes2(1, ruling_tris2(j, :)), ruling_nodes2(2, ruling_tris2(j, :)));
            if in == 1
                if on == 0
                    labels2(i, 1) = j;
                elseif isempty(ruling_idx)
                    labels2(i, 1) = j;
                else
                    labels2(i, 1) = nan;
                end
                break
            end
        end
    end
    labels2(labels2==0) = max(labels2(:, 1)); % bug workaround
    
    
    %% project unfolded onto a parabolic surface
    vert1_p = vert1;
    vert1_p(:, 3) = c.*(vert1_p(:, 1).^2 + vert1_p(:, 2).^2);
    
    vert2_p = vert2;
    vert2_p(:, 3) = c.*(vert2_p(:, 1).^2 + vert2_p(:, 2).^2);
    
    v_plane_p = v_plane;
    v_plane_p(3, :) = c.*(v_plane_p(1, :).^2 + v_plane_p(2, :).^2);
    
    m_plane_p = m_plane;
    m_plane_p(3, :) = c.*(m_plane_p(1, :).^2 + m_plane_p(2, :).^2);
    
    m_plane_ext_p = m_plane_ext;
    m_plane_ext_p(3, :) = c.*(m_plane_ext_p(1, :).^2 + m_plane_ext_p(2, :).^2);
    
    inner_edge_p = inner_edge;
    inner_edge_p(3, :) = c.*(inner_edge_p(1, :).^2 + inner_edge_p(2, :).^2);
    
    ruling_nodes1_p = [v_plane_p, m_plane_p(:, 2:end)];
    ruling_nodes2_p = [m_plane_ext_p, flip(inner_edge_p, 2), rot*v_plane_p(:, 1:end)];
    
    m_folded    = nan(3, length(v_plane_p));
    v_folded    = nan(3, length(v_plane_p));
    v_folded(:, 1) = [A; 0; min(v_plane_p(3, :))];
    m_folded(:, 1) = v_folded(:, 1);

    %% folded initial condition
    v_spiral = @(theta) [(A+(2*h)/beta*theta)*cos(theta); (A+(2*h)/beta*theta)*sin(theta); min(v_plane_p(3, :))];
    m_spiral = @(theta) [(A-h+(1-theta/beta)*h*(theta<beta)+(2*h)/beta*theta)*cos(theta); (A-h+(1-theta/beta)*h*(theta<beta)+(2*h)/beta*theta)*sin(theta); min(v_plane_p(3, :))];
    
    theta_v = 0;
    theta_m = 0;
    options = optimoptions('fmincon', 'Display', 'off');
    for i = 2:n
        dist_v = norm(v_plane_p(:, i-1)-v_plane_p(:, i));
        dist_m = norm(m_plane_p(:, i-1)-m_plane_p(:, i));
        fun_v = @(theta) (dist_v - norm(v_folded(:, i-1)-v_spiral(theta)))^2;
        th_v = fmincon(fun_v, theta_v(end)+beta, [], [], [], [], theta_v(end), theta_v(end)+beta, [], options);
        fun_m = @(theta) (dist_m - norm(m_folded(:, i-1)-m_spiral(theta)-[0;0;theta*delta_z/beta]))^2;
        th_m = fmincon(fun_m, theta_m(end)+beta, [], [], [], [], theta_m(end), theta_m(end)+beta, [], options);
        v_folded(:, i) = v_spiral(th_v);
        m_folded(:, i) = m_spiral(th_m) + [0;0;th_m*delta_z/beta];
        theta_v = [theta_v, th_v];
        theta_m = [theta_m, th_m];
    end
    
    m_folded_ext = m_folded;
    theta_ext_m = theta_m;
    
    for i = (n+1):n_ext
        dist_m = norm(m_plane_ext_p(:, i-1)-m_plane_ext_p(:, i));
        fun_m = @(theta) (dist_m - norm(m_folded_ext(:, i-1)-m_spiral(theta)-[0;0;theta*delta_z/beta]))^2;
        th_m = fmincon(fun_m, theta_ext_m(end)+beta, [], [], [], [], theta_ext_m(end), theta_ext_m(end)+beta, [], options);
        m_folded_ext(:, i) = m_spiral(th_m) + [0;0;th_m*delta_z/beta];
        theta_ext_m = [theta_ext_m, th_m];
    end
    
    z_v = zeros(size(theta_v));
    z_m = theta_m.*delta_z./beta;
    
    ruling_nodes1_f = [v_folded, m_folded(:, 2:end)];
    ruling_nodes2_f = [m_folded_ext, flip(inner_edge_p, 2), rot*v_folded(:, 1:end)];
    
    cone_idx = sum(theta_m < beta); % we won't constrain these nodes to the spiral
    
    ruling_nodes1_th = [theta_v, theta_m(2:end)];
    ruling_nodes1_z = [z_v, z_m(2:end)];
    
    ruling_nodes2_th = [theta_m, theta_v];
    ruling_nodes2_z = [z_m, z_v];
    
    %% Place subdivided nodes within ruling triangles
    
    vert1_f = zeros(length(vert1_p), 3);
    
    for i = 1:length(vert1_p)
        if ~isnan(labels1(i, 2)) 
            vert1_f(i, 1:3) = ruling_nodes1_f(:, labels1(i, 2));
        else
            tri_idx = labels1(i, 1);
            p = vert1_p(i, 1:2);
            p1 = [ruling_nodes1_p(1:2, ruling_tris1(tri_idx, 1)); 0];
            p2 = [ruling_nodes1_p(1:2, ruling_tris1(tri_idx, 2)); 0];
            p3 = [ruling_nodes1_p(1:2, ruling_tris1(tri_idx, 3)); 0];
            p1_f = ruling_nodes1_f(1:3, ruling_tris1(tri_idx, 1));
            p2_f = ruling_nodes1_f(1:3, ruling_tris1(tri_idx, 2));
            p3_f = ruling_nodes1_f(1:3, ruling_tris1(tri_idx, 3));
            
            dem = (p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*p3(1) + p2(2)*p3(1) + p1(1)*p3(2) - p2(1)*p3(2));
            beta_tri = -(p1(2)*p3(1) - p1(1)*p3(2) - p1(2)*p(1) + p3(2)*p(1) + p1(1)*p(2) - p3(1)*p(2))/dem;
            delta_tri = (p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*p(1) + p2(2)*p(1) + p1(1)*p(2) - p2(1)*p(2))/dem;
    
            vert1_f(i, 1:3) = p1_f + beta_tri*(p2_f-p1_f)+delta_tri*(p3_f-p1_f);
        end
    end
    
    vert2_f = zeros(length(vert2_p), 3);
    
    for i = 1:length(vert2_p)
        if ~isnan(labels2(i, 2)) 
            vert2_f(i, 1:3) = ruling_nodes2_f(:, labels2(i, 2));
        else
            tri_idx = labels2(i, 1);
            p = vert2(i, 1:2);
            p1 = [ruling_nodes2_p(1:2, ruling_tris2(tri_idx, 1)); 0];
            p2 = [ruling_nodes2_p(1:2, ruling_tris2(tri_idx, 2)); 0];
            p3 = [ruling_nodes2_p(1:2, ruling_tris2(tri_idx, 3)); 0];
            p1_f = ruling_nodes2_f(1:3, ruling_tris2(tri_idx, 1));
            p2_f = ruling_nodes2_f(1:3, ruling_tris2(tri_idx, 2));
            p3_f = ruling_nodes2_f(1:3, ruling_tris2(tri_idx, 3));
            
            dem = (p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*p3(1) + p2(2)*p3(1) + p1(1)*p3(2) - p2(1)*p3(2));
            beta_tri = -(p1(2)*p3(1) - p1(1)*p3(2) - p1(2)*p(1) + p3(2)*p(1) + p1(1)*p(2) - p3(1)*p(2))/dem;
            delta_tri = (p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*p(1) + p2(2)*p(1) + p1(1)*p(2) - p2(1)*p(2))/dem;
    
            vert2_f(i, 1:3) = p1_f + beta_tri*(p2_f-p1_f)+delta_tri*(p3_f-p1_f);
        end
    end
    
    %% combine gores into section
    
    vert_ref = vert1;
    vert_u = vert1_p;
    vert_f = vert1_f;
    tria_temp = zeros(size(tria2));
    labels = labels1(:, 2);
    
    for i = 1:length(vert2_p)
        found_idx = find(ismember(vert1_p, vert2_p(i, :), 'rows'));
        ruling_idx = find(ismembertol((rot*ruling_nodes1_p)', vert2_p(i, :), 1e-4, 'ByRows', true));
        tria_idx = find(tria2(:, 1:3)==i);
        if isempty(found_idx)
            vert_u = [vert_u; vert2_p(i, :)];
            vert_f = [vert_f; vert2_f(i, :)];
            vert_ref = [vert_ref; vert2(i, :)];
            tria_temp(tria_idx) = length(vert_u);
            if ~isempty(ruling_idx)
                labels = [labels; -ruling_idx];
            else
                labels = [labels; nan];
            end
        else
            tria_temp(tria_idx) = found_idx;
        end
        
    end
    
    faces = [tria1; tria_temp];
    
    % assign edges
    for i = 1:length(faces)
        edge1 = faces(i, 1:2);
        edge2 = faces(i, 2:3);
        edge3 = faces(i, [1 3]);
        if i == 1
            edges = [edge1; edge2; edge3];
        else
            if isempty(find(ismember(edges, edge1, 'rows'))) & isempty(find(ismember(edges, flip(edge1), 'rows')))
                edges = [edges; edge1];
            end
            if isempty(find(ismember(edges, edge2, 'rows'))) & isempty(find(ismember(edges, flip(edge2), 'rows')))
                edges = [edges; edge2];
            end
            if isempty(find(ismember(edges, edge3, 'rows'))) & isempty(find(ismember(edges, flip(edge3), 'rows')))
                edges = [edges; edge3];
            end
        end
    end

    %% find indices of constrained nodes

    inner_idx = find(ismember(vert_f(:, 1:3), v_folded(1:3, 1)', 'rows'));
    inner_idx_rot = find(ismember(vert_f(:, 1:3), (rot*v_folded(1:3, 1))', 'rows'));
    stat_idxs = [inner_idx, inner_idx_rot];
    vert_f(stat_idxs(end), :) = vert_u(stat_idxs(end), :);
    for i = 1:length(div3)
        stat_idx = find(ismember(vert_f(:, 1:3), (inner_edge_p(1:3, i))', 'rows'));
        stat_idxs = [stat_idxs, stat_idx];
        vert_f(stat_idxs(end), :) = vert_u(stat_idxs(end), :);
    end
    
    outer_idx = [];
    outer_idx_rot = find(ismember(vert_f(:, 1:3), (rot*v_folded(1:3, end))', 'rows'));
    for i = 1:length(vert_u)
        vert_temp = vert_u(i, 1:2);
        if ismembertol(norm(vert_temp), R, 1e-7) && i ~= outer_idx_rot
            outer_idx = [outer_idx, i];
        end
    end
    % outer_idx1 = find(ismember(vert_f(:, 1:3), v_folded(1:3, end)', 'rows'));
    % outer_idx2 = find(ismember(vert_f(:, 1:3), m_folded(1:3, end)', 'rows'));
    % 
    % outer_idxs = 
    % for i = 1:length(end_n-2)

        %% ribs (if non-zero rib depth specified)
    if rib_d > 0
        nNodes = length(vert_u);
        % unfolded
        rib_n = max(round(rib_d/((R-A)/n))+1, 2);
        rib_nodes_u = zeros(3, length(v_plane_p)*(rib_n-1)); 
        rib_z = linspace(0, rib_d, rib_n);
        rib_edges = [];
        rib_faces = [];
        for j = 1:(rib_n-1)
            rib_nodes_u(:, ((j-1)*length(v_plane_p)+1):(j*length(v_plane_p))) = v_plane_p - [zeros(2, length(v_plane_p)); rib_z(j+1).*ones(1, length(v_plane_p))];
            
        end

        % folded
        rib_nodes_f = zeros(3, length(v_folded)*(rib_n-1)); 
        for j = 1:(rib_n-1)
            rib_nodes_f(:, ((j-1)*length(v_folded)+1):(j*length(v_folded))) = v_folded - [zeros(2, length(v_folded)); rib_z(j+1).*ones(1, length(v_folded))];
        end

        % ref
        rib_nodes_ref = zeros(3, length(v_plane)*(rib_n-1)); 
        for j = 1:(rib_n-1)
            rib_nodes_ref(:, ((j-1)*length(v_plane)+1):(j*length(v_plane))) = v_plane - [zeros(2, length(v_plane)); rib_z(j+1).*ones(1, length(v_plane))];
        end
        
        rib_edges = [];
        rib_faces = [];
        for k = 1:(n-1) % first row
            idx1 = find(labels == k);
            idx2 = find(labels == k+1);
            rib_edges = [rib_edges; idx1, nNodes+k; idx1, nNodes+k+1];
            rib_faces = [rib_faces; idx1, idx2, nNodes+k+1; idx1, nNodes+k+1, nNodes+k];
        end
        rib_edges = [rib_edges; find(labels == n), nNodes+n];
        for j = 1:(rib_n-2)
            for k = 1:(n-1)
                base_idx = nNodes + (j-1)*n;
                rib_edges = [rib_edges; base_idx+k, base_idx+k+1; base_idx+k, base_idx+k+n; base_idx+k, base_idx+k+n+1];
                rib_faces = [rib_faces; base_idx+k, base_idx+k+1, base_idx+k+n+1; base_idx+k, base_idx+k+n+1, base_idx+k+n];
            end
            rib_edges = [rib_edges; nNodes+j*n, nNodes+(j+1)*n]; % end edge
        end
        for k = 1:(n-1)
            rib_edges = [rib_edges; nNodes+(rib_n-2)*n+k, nNodes+(rib_n-2)*n+k+1]; % botton edge
        end

        vert_u = [vert_u; rib_nodes_u'];
        vert_f = [vert_f; rib_nodes_f'];
        vert_ref = [vert_ref; rib_nodes_ref'];
        edges = [edges; rib_edges];
        faces = [faces; rib_faces];
        labels = [labels; nan(length(rib_nodes_u), 1)];
        stat_idxs = [stat_idxs, nNodes+1+n.*(0:(rib_n-2))];
        outer_idx = [outer_idx, nNodes+n+n.*(0:(rib_n-2))];
    end
%%
    % fig = figure;
    % patch('faces',faces(:,1:3),'vertices',vert_u(:, 1:3), ...
    %         'facecolor',[0.5, 0.5, 0.5], 'facealpha', 0.5, ...
    %         'edgecolor',[1,0,0], 'edgealpha', 0.2) ;
    % hold on
    % plot3dNodesEdges(vert_u', edges, nan(length(edges), 1), fig);


    %% Plot

    if plot_on
        % plot unfolded
        
        figure('Color', [1 1 1]);
        inner = [];
        patch('faces',faces(:,1:3),'vertices',vert_u(:, 1:2), ...
            'facecolor','w', 'facealpha', 0, ...
            'edgecolor',[1,0,0], 'edgealpha', 0.2) ;
        hold on; axis image off;
        
        patch('faces',ruling_tris1(:,1:2),'vertices', ruling_nodes1(1:2, :)', ...
            'facecolor','w', ...
            'edgecolor',[.1,.1,.1], ...
            'linewidth',1.5) ;
        gscatter(vert1_p(:, 1), vert1_p(:, 2), labels1(:, 1))
        
        % patch('faces',tria2(:,1:3),'vertices',vert2(:, 1:2), ...
        %     'facecolor','w', 'facealpha', 0, ...
        %     'edgecolor',[0,0,1], 'edgealpha', 0.2) ;
        patch('faces',ruling_tris2(:,1:2),'vertices', ruling_nodes2(1:2, :)', ...
            'facecolor','w', ...
            'edgecolor',[.1,.1,.1], ...
            'linewidth',1.5) ;
        gscatter(vert2(:, 1), vert2(:, 2), labels2(:, 1))
        
        v_verts = v_plane(1:2, :);
        m_verts = m_plane(1:2, :);
        plot(v_verts(1, :), v_verts(2, :), 'b-', 'LineWidth', 1.5)
        plot(v_verts(1, :), v_verts(2, :), 'k.', 'MarkerSize', 1)
        plot(m_verts(1, :), m_verts(2, :), 'r-', 'LineWidth', 1.5)
        plot(m_verts(1, :), m_verts(2, :), 'k.', 'MarkerSize', 1)
        inner = [inner, v_plane(1:2, 1)];
        
        inner = [inner, inner(:, 1)];
        plot(inner(1, :), inner(2, :), 'k-', 'LineWidth', 1.5)
        axis equal
        axis tight
        axis off
        
        ax = gca;
        ax.Clipping = 'off';
        
        % plot folded
        figure('Color', [1 1 1]);
        % plot3(v_folded(1, :), v_folded(2, :), v_folded(3, :), 'b-')
        % hold on
        % plot3(m_folded(1, :), m_folded(2, :), m_folded(3, :), 'r-')
        
        inner = [];
        patch('faces',faces(:,1:3),'vertices',vert_f(:, 1:3), ...
            'facecolor','w', 'facealpha', 0, ...
            'edgecolor',[1,0,0], 'edgealpha', 0.2) ;
        patch('faces',ruling_tris2(:,1:2),'vertices', ruling_nodes2_f(1:3, :)', ...
            'facecolor','w', ...
            'edgecolor',[.1,.1,.1], ...
            'linewidth',1.5) ;
        hold on; axis image off;
        
        % patch('faces',ruling_tris1(:,1:3),'vertices',ruling_nodes1_f(1:3, :)', ...
        %     'facecolor','w', ...
        %     'edgecolor',[.1,.1,.1], ...
        %     'linewidth',1.5) ;
        % gscatter(vert1_f(:, 1), vert1_f(:, 2), vert1_f(:, 3))
        % 
        % patch('faces',tria2(:,1:3),'vertices',vert2_f(:, 1:3), ...
        %     'facecolor',[0.7 0.7 0.7], 'facealpha', 0, ...
        %     'edgecolor',[0,0,1], 'edgealpha', 0.3) ;
        
        
        v_verts = v_folded(1:3, :);
        m_verts = m_folded(1:3, :);
        plot3(v_verts(1, :), v_verts(2, :), v_verts(3, :), 'b-', 'LineWidth', 1.5)
        plot3(v_verts(1, :), v_verts(2, :), v_verts(3, :), 'k.', 'MarkerSize', 1)
        plot3(m_verts(1, :), m_verts(2, :), m_verts(3, :), 'r-', 'LineWidth', 1.5)
        plot3(m_verts(1, :), m_verts(2, :), m_verts(3, :), 'k.', 'MarkerSize', 1)
        % inner = [inner, v_plane(1:2, 1)];
        % 
        % inner = [inner, inner(:, 1)];
        % plot(inner(1, :), inner(2, :), 'k-', 'LineWidth', 1.5)
        axis equal
        axis tight
        axis off
        
        ax = gca;
        ax.Clipping = 'off';
    end

end

function alphas = getAlphas(nodes, faces)
    nodes = nodes;
    alphas = zeros(size(faces));
    for i = 1:length(faces)
        p1_index = faces(i, 1);
        p2_index = faces(i, 2);
        p3_index = faces(i, 3);

        p1 = nodes(p1_index, :);
        p2 = nodes(p2_index, :);
        p3 = nodes(p3_index, :);

        alpha1 = angleBetweenVectors3d(p2-p1, p3-p1);
        alpha2 = angleBetweenVectors3d(p1-p2, p3-p2);
        alpha3 = angleBetweenVectors3d(p1-p3, p2-p3);

        alphas(i, :) = [alpha1 alpha2 alpha3];
    end
end