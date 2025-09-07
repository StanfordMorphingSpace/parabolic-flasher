clc; clear;

tic

% add all subfolders
addpath(genpath(pwd))

save_on = 1;
plot_on = 1;
save_path = "D:\Curved_crease_antennas\SciTech_2026\ribs_results";

% wildtronics dish
A = 115/2/1000; % 4 in
N = 8;
h = 0.8/1000; % layer thickness
n = 60; % total subdivisions
R = 507/2/1000; % outer radius as measured
c = 0.25*1/(4*118.11/1000); %6e-2 % 4.65in focus to vertex
%c = 0;
iter = 10000;
rib_d = 0.015;

l = (R-A)/n/sqrt(3);
rib_n = max(round(rib_d/((R-A)/n))+1, 2);

% % FLUTE
% A = 1.5; % 
% N = 18;
% h = 0.015; % layer thickness
% n = 50; % total subdivisions
% R = 50/2; % outer radius
% c = 1/(4*(2*R*1.5)); %6e-2 % 4.65in focus to vertex
% iter = 10000;
% l = (R-A)/n/sqrt(3);

% layer height 0.2
% A = 20; % 
% N = 6;
% h = 0.8; % layer thickness
% n = 24; % total subdivisions
% R = 120; % 
% c = 4e-3; % 4.65in focus to vertex

geo.A = A;
geo.N = N;
geo.h = h; % layer thickness
geo.n = n; % total subdivisions
geo.R = R; % outer radius
geo.c = c;

% % material properties polycarbonate
E = 2390000000; 
v = 0.37;
t = 0.0006; % 0.00076 reported, 0.0006 measured
% %t = 0.001;

% material properties carbon fiber
% E = 120E9; 
% v = 0.28;
% t = 0.01; % 

geo.k_axial = 8*E*t*l/(1-v^2)/9;

geo.k_fold  = E*t^3/(12*l*(1-v^2)); % 0.7
geo.k_face  = geo.k_fold; % 0.2
% k_facet = 0.7;
geo.gamma   = 0.3; % damping 0.45 0.2

dt = 0.5*1/(2*pi*sqrt(geo.k_axial/l*(sqrt(3)/2)));
dt = 0.5*1/(2*pi*sqrt(geo.k_axial));
mass_scalar = 1e3;

% geo.k_axial = 40;
% geo.k_face  = 0.2; % 0.2
% geo.k_fold  = 0.1; % 0.7
% % k_facet = 0.7;
% geo.gamma   = 0.2; % damping 0.45 0.2

beta    = 2*pi()/N;
rot     = [ cos(beta), -sin(beta), 0;...
            sin(beta), cos(beta), 0;...
            0, 0, 1];

[vert_u, vert_f, vert_ref, labels, edges, faces, stat_idxs, outer_idx, cone_idx] = generateMesh_ribs(A, N, h, n, R, c, 0, rib_d);
%%
adj_face = cell(length(edges), 1);
% rearrange how faces are stored
for i = 1:length(edges)
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);
        
        % Bending forces
        faces_with_node1        = (faces(:, 1) == p1_index) | (faces(:, 2) == p1_index) | (faces(:, 3) == p1_index);
        faces_with_node2        = (faces(:, 1) == p2_index) | (faces(:, 2) == p2_index) | (faces(:, 3) == p2_index);
        faces_adjacent_to_edge  = faces((faces_with_node1 & faces_with_node2), :);
        
        adj_face{i} = faces_adjacent_to_edge(~ismember(faces_adjacent_to_edge,[p1_index p2_index]));
end

% Only keep edges with 2 adjacent faces
mask = cellfun(@(x) numel(x)==2, adj_face);
creaseEdges = edges(mask,:);
creaseIdx   = find(mask);

if isempty(creaseIdx)
    F_crease = zeros(nNodes,3);
    E_crease = zeros(nEdges,1);
    return;
end

% Build adjacency matrix safely
adj_mat = cell2mat(cellfun(@(x) x(:).', adj_face(mask), 'UniformOutput', false));

adj_faces.creaseIdx = creaseIdx;
adj_faces.adj_mat = adj_mat;
adj_faces.mask = mask;

nNodes = length(vert_u);

% fix normals
focal_point = [0, 0, 1/(c*4)];

p1    = vert_u(creaseEdges(:,1), :);
p2    = vert_u(creaseEdges(:,2), :);
p3_1  = vert_u(adj_mat(:,1), :);
p3_2  = vert_u(adj_mat(:,2), :);

e = p2 - p1;
e = e ./ vecnorm(e,2,2);

n1 = cross(p2 - p1, p3_1 - p1, 2);
n2 = cross(p3_2 - p1, p2 - p1, 2);

n1 = n1 ./ vecnorm(n1,2,2);
n2 = n2 ./ vecnorm(n2,2,2);

ref = focal_point-p1;
creaseEdges(dot(ref, n1, 2) < 0, :) = flip(creaseEdges(dot(ref, n1, 2) < 0, :), 2);
adj_faces.creaseEdges = creaseEdges;

%%
lengths     = getEdgeLengths(vert_u(:, 1:3), edges);
lengths_p   = getEdgeLengths(vert_f, edges);
error       = (lengths - lengths_p)./lengths_p*100;

figure();
plot(1:length(edges), error, 'o--')
xlabel('Edge index')
ylabel('Length error (percent)')

angles_f = foldedCreaseAngles_fast(vert_f, vert_u, edges, adj_faces);
angles_u = foldedCreaseAngles_fast(vert_u, vert_f, edges, adj_faces);
angles = mod((angles_f - angles_u) + pi, 2*pi) - pi;
angles = sign(angles)*pi - angles;
angles(angles==0) = pi;

fig_deployed = plot3dNodesEdges(vert_u', edges, angles); %nan(length(edges), 1)
figure(fig_deployed)

%% Dynamic Relaxation

mass = mass_scalar*ones(length(vert_u), 1);

k_fold = geo.k_fold*ones(length(edges), 1);
% for i = 1:length(edges)
%     p1 = edges(i, 1);
%     p2 = edges(i, 2);
%     if labels(p1) > n && labels(p2) > n % mountain
%         k_fold(i) = 0;
%     elseif (labels(p1) == 1 && labels(p2) > n) || (labels(p2) == 1 && labels(p1) > n)
%         k_fold(i) = 0;
%     end
% end

for i = 1:length(edges)
    p1 = edges(i, 1);
    p2 = edges(i, 2);
    if ~isnan(labels(p1)) && ~isnan(labels(p2))
        k_fold(i) = 0;
    end
end

k_fold(~mask) = 0;

geo.k_fold = k_fold;

p_f = vert_f(:, 1:3);
v_f = zeros(nNodes, 3); % initial velocity = 0

p_u = vert_u(:, 1:3);
v_u = zeros(nNodes, 3); % initial velocity = 0

vert_u = zeros(nNodes, 3, iter+1);
vert_u(:, :, 1) = p_u;
vert_f = zeros(nNodes, 3, iter+1);
vert_f(:, :, 1) = p_f;

accel_u = zeros(nNodes, 3, iter+1);
accel_f = zeros(nNodes, 3, iter+1);

E_cr = zeros(1, iter);
E_ax = zeros(1, iter);
E_v = zeros(1, iter);
E_tot_prev = Inf;
E_diff = Inf;
for i = 1:iter
    if E_diff < 1e-9
        disp("converged")
        break
    end
    [p_u, p_f, v_u, v_f, a_u, a_f, E_crease, E_axial, E_v_u, E_v_f] = DR_Step_parabolic(p_u, p_f, v_u, v_f, vert_ref, labels, edges, adj_faces, dt, geo, mass, i, stat_idxs, outer_idx, R, cone_idx, n, rib_n, rib_d);
    
    % Store Energy
    E_cr(i) = sum(E_crease);
    E_ax(i) = sum(E_axial);
    E_v(i) = sum(E_v_f) + sum(E_v_u);
    
    E_diff = abs(E_cr(i)+E_ax(i)+E_v(i)-E_tot_prev);
    E_tot_prev = E_cr(i)+E_ax(i)+E_v(i);

    vert_f(:, :, i+1) = p_f;
    accel_f(:, :, i+1) = a_f;

    vert_u(:, :, i+1) = p_u;
    accel_u(:, :, i+1) = a_u;
end

end_incr = i;

toc

%% Plot Energies

if plot_on

    figure;
    plot(1:length(E_ax), E_ax, "LineWidth", 2)
    set(gca, 'YScale', 'log')
    xlabel("Iterations")
    ylabel("Total Stretching Energy [J]")
    grid on
    set(gca, "FontSize", 18)
    
    figure;
    plot(1:length(E_cr), E_cr, "LineWidth", 2)
    set(gca, 'YScale', 'log')
    xlabel("Iterations")
    ylabel("Total Bending Energy [J]")
    grid on
    set(gca, "FontSize", 18)
    
    
    figure;
    plot(1:length(E_v), E_v, "LineWidth", 2)
    set(gca, 'YScale', 'log')
    xlabel("Iterations")
    ylabel("Total Kinetic Energy [J]")
    grid on
    set(gca, "FontSize", 18)
    
    
    figure;
    plot(1:length(E_v), E_ax+E_cr+E_v, "LineWidth", 2)
    set(gca, 'YScale', 'log')
    xlabel("Iterations")
    ylabel("Total Energy [J]")
    grid on
    set(gca, "FontSize", 18)
end


%% Get edge lengths and angles
lengths     = getEdgeLengths(vert_u(:, 1:3, end), edges);
lengths_p   = getEdgeLengths(vert_f(:, :, end), edges);
error       = (lengths - lengths_p)./lengths_p*100;

figure();
plot(1:length(edges), error, 'o--')
xlabel('Edge index')
ylabel('Length error (percent)')

angles_f = foldedCreaseAngles_fast(vert_f(:, :, end), vert_u(:, :, end), edges, adj_faces);
angles_u = foldedCreaseAngles_fast(vert_u(:, :, end), vert_f(:, :, end), edges, adj_faces);
angles = mod((angles_f - angles_u) + pi, 2*pi) - pi;
angles = sign(angles)*pi - angles;
angles(angles==0) = pi;

%% plotting

if plot_on
    % plot deployed
    deployed = figure('Color', [1 1 1]);
    hold on
    for i = 0:(n-1)
        %deployed = plot3dNodesEdges((rot^i*vert_u(:, :, end)'), edges, angles, deployed);
        patch('faces',faces(:,1:3),'vertices',(rot^i*vert_u(:, 1:3, end)')', ...
            'facecolor',[0.7 0.7 0.7], 'facealpha', 0.6, ...
            'edgecolor',[1,0,0], 'edgealpha', 1) ;
    end
    hold off
    axis equal; axis tight; axis off
    
    ax = gca; ax.Clipping = 'off';
    
    % plot folded
    stowed = figure('Color', [1 1 1]);
    inner = [];
    for i = 0:(n-1)
        %stowed = plot3dNodesEdges((rot^i*vert_f(:, 1:3, end)'), edges, angles, stowed);
        hold on
        patch('faces',faces(:,1:3),'vertices',(rot^i*vert_f(:, 1:3, end)')', ...
            'facecolor',[0.7 0.7 0.7], 'facealpha', 0.4, ...
            'edgecolor',[1,0,0], 'edgealpha', 1) ;
    end
    hold off
    axis equal; axis tight; axis off
    
    ax = gca; ax.Clipping = 'off';
end

%% Save major fold lines

if save_on
    major_v_u = vert_u(0 < labels & labels <= n, :, end);
    order = labels(0 < labels & labels <= n)-min(labels(0 < labels & labels <= n))+1;
    major_v_u = flip(sortrows([major_v_u, order], 4), 1);
    %major_v_u = flip(major_v_u(order, :))
    major_v_f = vert_f(labels(0 < labels & labels <= n), :, end);
    major_v_f = flip(sortrows([major_v_f, order], 4), 1);
    %major_v_f = flip(major_v_f(order, :))
    
    major_m_u = vert_u((labels > n), :, end);
    order = labels(labels > n)-min(labels(labels>n))+1;
    major_m_u = sortrows([major_m_u, order], 4);
    %major_m_u = major_m_u(order, :)
    major_m_f = vert_f((labels > n), :, end);
    major_m_f = sortrows([major_m_f, order], 4);
    %major_m_f = major_m_f(order, :)
    
    save(fullfile(save_path, sprintf("major_folds_quarter_n%d_N%d_rib%d.mat", [n, N, rib_d*1000])), 'major_v_u', 'major_v_f', 'major_m_u', 'major_m_f');
    writematrix(major_v_u(:, 1:2), fullfile(save_path, sprintf('major_v_u_quarter_n%d_N%d_rib%d.csv', [n, N, rib_d*1000])));
    writematrix(major_m_u(:, 1:2), fullfile(save_path, sprintf('major_m_u_quarter_n%d_N%d_rib%d.csv', [n, N, rib_d*1000])));
    
    inner_idxs = vert_u(stat_idxs([1, end:-1:2]), :, end);
    hexagon = inner_idxs;
    
    for i = 1:(N-1)
        hexagon = [hexagon; (rot^i*inner_idxs')'];
    end
    
    hexagon = [hexagon; hexagon(1, :)];
    
    writematrix(hexagon(:, 1:2), fullfile(save_path, sprintf('inner_hexagon_quarter_n%d_N%d_rib%d.csv', [n, N, rib_d*1000])));
    save(fullfile(save_path, sprintf("120524_converge_quarter_n%d_N%d_rib%d.mat", [n, N, rib_d*1000])));
end

%% Generate unified mesh (This takes a while)
% [unfolded, edges_one, faces_one] = makeFullMesh(vert_u(:, :, end), edges, faces, rot, N);
% [folded, ~, ~] = makeFullMesh(vert_f(:, :, end), edges, faces, rot, N);
% 
% angles_one = foldedCreaseAngles(unfolded, folded, edges_one, faces_one);
% angles_one = sign(angles_one)*pi - angles_one;
% angles_one(angles_one==0) = pi;
%
% %%
% folded_one = figure('Color', [1 1 1]);
% patch('faces',faces_one,'vertices',folded, ...
%         'facecolor',[0.7 0.7 0.7], 'facealpha', 0.4, ...
%         'edgecolor',[1,0,0], 'edgealpha', 0) ;
% 
% folded_one = plot3dNodesEdges(folded', edges_one, angles_one, folded_one);
% 
% axis equal; axis tight; axis off
% 
% ax = gca; ax.Clipping = 'off';
% 
% deployed_one = figure('Color', [1 1 1]);
% patch('faces',faces_one,'vertices',unfolded, ...
%         'facecolor',[0.7 0.7 0.7], 'facealpha', 0.5, ...
%         'edgecolor',[1,0,0], 'edgealpha', 0) ;
% 
% deployed_one = plot3dNodesEdges(unfolded', edges_one, angles_one, deployed_one);
% 
% axis equal; axis tight; axis off; view(-37.5, 30)
% 
% ax = gca; ax.Clipping = 'off';
% 
% if save_on
%     triag = triangulation(faces_one, unfolded);
%     stlwrite(triag, fullfile(save_path, sprintf('DR_parabolic_test_v4_N8_n%d_c0.stl', n)));
% end
%% plot folded with time slider (Out of date)

% S.vert_f = vert_f;
% S.vert_u = vert_u;
% S.faces = faces;
% S.ind = end_incr;
% S.accel_f = accel_f;
% S.accel_u = accel_u;
% S.rot = rot;
% 
% S.ruling_nodes1_f = ruling_nodes1_f;
% S.n = n;
% 
% S.fig = figure('Units', 'normalized', 'Position', [0.05 0.01 0.9 0.9]);
% xlim([-0.5, 1.5]);
% ylim([-0.5, 1.5]);
% zlim([-0.5, 1.5]);
% plot_ind(S);
% 
% 
% uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.04, 0.01, 0.05, 0.04], 'String', 'ind')
% shiftXSlider = uicontrol(   'Style', 'slider',...
%                             'Units', 'normalized', 'Position', [0.1, 0.03, 0.8, 0.03],...
%                             'value', S.ind, 'min', 1, 'max', end_incr,...
%                             'callback', {@showEdgeworkSliderCB, 'ind'});
% guidata(S.fig, S)

%% slider functions

function showEdgeworkSliderCB(slider, EventData, Param)    
    S = guidata(slider);
    S.(Param) = get(slider, 'Value');
    S = plot_ind(S);
    guidata(slider, S);
end

function S = plot_ind(S)
    vert_f = S.vert_f(:, :, floor(S.ind));
    accel_f = S.accel_f(:, :, floor(S.ind));
    vert_u = S.vert_u(:, :, floor(S.ind));
    accel_u = S.accel_u(:, :, floor(S.ind));
    S.xlim = xlim;
    S.ylim = ylim;
    S.zlim = zlim;
    cla;
    
    patch('faces',S.faces(:,1:3),'vertices',vert_f(:, 1:3), ...
        'facecolor','w', 'facealpha', 0, ...
        'edgecolor',[1,0,0], 'edgealpha', 0.4) ;
    hold on; axis image off;
    patch('faces',S.faces(:,1:3),'vertices',vert_u(:, 1:3), ...
        'facecolor','w', 'facealpha', 0, ...
        'edgecolor',[0,0,1], 'edgealpha', 0.4) ;
    
    v_verts = S.ruling_nodes1_f(1:3, 1:S.n);
    v_verts_rot = S.rot*S.ruling_nodes1_f(1:3, 1:S.n);
    m_verts = S.ruling_nodes1_f(1:3, (S.n+1):end);
    plot3(v_verts(1, :), v_verts(2, :), v_verts(3, :), 'b-', 'LineWidth', 1.5)
    plot3(v_verts(1, :), v_verts(2, :), v_verts(3, :), 'k.', 'MarkerSize', 1)
    plot3(m_verts(1, :), m_verts(2, :), m_verts(3, :), 'r-', 'LineWidth', 1.5)
    plot3(m_verts(1, :), m_verts(2, :), m_verts(3, :), 'k.', 'MarkerSize', 1)
    plot3(v_verts_rot(1, :), v_verts_rot(2, :), v_verts_rot(3, :), 'b-', 'LineWidth', 1.5)
    plot3(v_verts_rot(1, :), v_verts_rot(2, :), v_verts_rot(3, :), 'k.', 'MarkerSize', 1)

    quiver3(vert_f(:, 1),vert_f(:, 2),vert_f(:, 3), accel_f(:, 1),accel_f(:, 2),accel_f(:, 3))
    quiver3(vert_u(:, 1),vert_u(:, 2),vert_u(:, 3), accel_u(:, 1),accel_u(:, 2),accel_u(:, 3))
    axis equal
    axis tight
    axis off
    xlim(S.xlim)
    ylim(S.ylim)
    zlim(S.zlim)
    
    ax = gca;
    ax.Clipping = 'off';
end

