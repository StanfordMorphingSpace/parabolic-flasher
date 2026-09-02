clc; clear;

% add all subfolders
addpath(genpath(pwd))

save_on = 0; % toggle save fold lines for abaqus
plot_on = 1; % toggle plotting figures
save_stability = 0; % toggle saving convergence info
save_path = "D:\Curved_crease_antennas\SciTech_2027\fold_pattern\brims";

%% wildtronics dish
A = 115/2/1000; % inner polygon radium, m (4 in)
N = 8;
h = 2.2/1000; % layer thickness, m
n = 50; % total subdivisions per major fold line
R = 507/2/1000; % outer radius as measured, m
c = 1/(4*118.11/1000); % 4.65in focus to vertex for test article, m
c = 0.3*c; % scale c for different depths
iter = 5; % number of fmincon iterations

surf_func = @(r) c*r.^2; % surface function for the paraboloid
surf_func_prime = @(r) 2*c*r; % d(surf_func)/dr

% no rib set rib_d = 0
rib_d = 0.00; % dib depth, m
% no brim set brim_R = 0
brim_R = 0.554/2; % brim outer radius

brim_R = 0;

brim_a = (c*R*(R-2*brim_R) + 0.142)/(brim_R-R)^2;
brim_b = 2*R*(c*brim_R^2 - 0.142)/(brim_R-R)^2;
brim_c = R^2*(0.142 - c*brim_R^2)/(brim_R-R)^2;

brim_func = @(r) brim_a*r.^2 + brim_b*r + brim_c;
brim_func_prime = @(r) 2*brim_a*r + brim_b;


% % material properties polycarbonate
E = 2390000000; 
v = 0.37;
t = 0.002; % 0.00076 reported, 0.0006 measured outside, 0.0007 measured inside

% A = 0.025; % inner polygon radium, m (4 in)
% N = 10;
% h = 0.6/1000; % layer thickness, m
% n = 100; % total subdivisions per major fold line
% R = 0.5; % outer radius as measured, m
% c = 1/(4*1); % 4.65in focus to vertex for test article, m
% %c = 0.*c; % scale c for different depths
% iter = 50000; % number of DR iterations (~30000 for convergence)
% rib_d = 0.00; % dib depth, m
% 
% % % material properties polycarbonate
% E = 2390000000; 
% v = 0.37;
% t = 0.0005; % 0.00076 reported, 0.0006 measured

% % %% Lunar Reflector prototype
% A = 0.08/2; % inner polygon radius, m
% N = 8; % rotational symmetry
% h = 1.2*0.67/1000; % layer thickness, m
% n = 55; % total subdivisions per major fold line
% R = 0.4318/2; % outer radius as measured, m
% c = 1/(4*1); % 1m focal length, m
% iter = 100000; % number of DR iterations (~30000 for convergence)
% rib_d = 0; % dib depth, m
% 
% % 
% % % material properties carbon fiber
% E = 120E9; 
% v = 0.28;
% t = 0.67/1000; % material thickness
% %t = 2*2.54e-4;

% %% Parachute
% A = 1.505/2; % inner polygon radium, m (4 in)
% N = 8; % parachute has 80 gores
% h = 0.002; % layer separation (should be greater than t), m
% n = 60; % total subdivisions per major fold line
% R = 15.652/2; % outer radius as measured, m
% f_num = 1; % f/D ratio
% c = 1/(4*f_num*2*R); % chord length, m
% c = 0*c; % scale c for different depths
% iter = 30000; % number of DR iterations (~30000 for convergence)
% rib_d = 0; % dib depth, m
% 
% % material properties parachute
% E = 9e9; 
% v = 0.4;
% t = 0.000076; % 0.00076 reported, 0.0006 measured

rib_n = max(round(rib_d/((R-A)/n))+1, 2);
if rib_d == 0
    rib_n = 0;
end

brim_n = max(round((brim_R-R)/((R-A)/n))+1, 2);
if brim_R == 0
    brim_n = 0;
end

% % FLUTE
% A = 1.5; % 
% N = 18;
% h = 0.015; % layer thickness
% n = 50; % total subdivisions
% R = 50/2; % outer radius
% c = 1/(4*(2*R*1.5)); %6e-2 % 4.65in focus to vertex
% iter = 10000;
% l = (R-A)/n/sqrt(3);

geo.A = A;
geo.N = N;
geo.h = h; % layer thickness
geo.n = n; % total subdivisions
geo.R = R; % outer radius
geo.surf_func = surf_func;
geo.surf_func_prime = surf_func_prime;
geo.brim_R = brim_R; % brim outer radius
geo.brim_func = brim_func;
geo.brim_func_prime = brim_func_prime;
geo.E = E;
geo.v = v;
geo.rib_d = rib_d;
geo.t = t;

beta    = 2*pi()/N;
rot     = [ cos(beta), -sin(beta), 0;...
            sin(beta), cos(beta), 0;...
            0, 0, 1];

% generate mesh
%[vert_u, vert_f, vert_ref, labels, edges, faces] = generateMesh_ribs(A, N, h, n, R, surf_func, 0, rib_d, brim_R, brim_func);
[vert_u, vert_f, vert_ref, labels, edges, faces] = generateMesh_ribs_old(A, N, h, n, R, surf_func, 0, rib_d);


lengths = getEdgeLengths(vert_u, edges);

l = mean(lengths)*sqrt(3)/2;

%l2 = (R-A)/n/sqrt(3);
geo.k_axial = 8*E*t*l/(1-v^2)/9;

geo.k_fold  = E*t^3/(12*l*(1-v^2)); % 0.7
geo.gamma   = 0.1; % damping 0.45 0.2

mass_scalar = 1;

%dt = 1/(2*pi*sqrt(geo.k_axial/(l*(2/sqrt(3)))/mass_scalar));
dt = 0.9*1/(2*pi*sqrt(geo.k_axial/min(lengths)/mass_scalar));
%dt = 0.5*1/(2*pi*sqrt(geo.k_axial));

%%
adj_faces = getAdjacentFaces(vert_u, edges, faces);

nNodes = length(vert_u);

%% Initial error and angles
lengths     = getEdgeLengths(vert_u(:, 1:3), edges);
lengths_p   = getEdgeLengths(vert_f, edges);
error       = (lengths - lengths_p)./lengths*100;

figure();
plot(1:length(edges), error, 'o--')
xlabel('Edge index')
ylabel('Length error (percent)')

angles_f = foldedCreaseAngles_fast(vert_f, vert_u, edges, adj_faces);
angles_u = foldedCreaseAngles_fast(vert_u, vert_f, edges, adj_faces);
angles = angles_f - angles_u + pi;

fig_deployed = plot3dNodesEdges(vert_u', edges, angles); %nan(length(edges), 1)
figure(fig_deployed)

%% Dynamic Relaxation
tic

k_fold = geo.k_fold*ones(length(edges), 1);

% zero out fold stiffness along ruling-line edges (they aren't actual creases)
ruled_idx = [labels.valley_idx, labels.valley_rot_idx, labels.mount_idx];
if isfield(labels, 'ribs_idx')
    ruled_idx = [ruled_idx, labels.ribs_idx];
end
if isfield(labels, 'brim_idx')
    ruled_idx = [ruled_idx, labels.brim_valley_idx, labels.brim_valley_rot_idx, labels.brim_mountain_idx];
end

for i = 1:length(edges)
    p1 = edges(i, 1);
    p2 = edges(i, 2);
    if ismember(p1, ruled_idx) && ismember(p2, ruled_idx)
        k_fold(i) = 0;
    elseif ismember(p1, labels.outer_idx) && ismember(p2, labels.outer_idx)
        k_fold(i) = 0;
    end

end

k_fold(~adj_faces.mask) = 0;

geo.k_fold = k_fold;

% unity mass matrix
mass = mass_scalar*ones(length(vert_u), 1);

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
    if E_diff < 5e-7
        disp("converged")
        break
    end
    [p_u, p_f, v_u, v_f, a_u, a_f, E_crease, E_axial, E_v_u, E_v_f] = DR_Step_parabolic(p_u, p_f, v_u, v_f, labels, edges, adj_faces, dt, geo, mass, i);
    
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

runtime = toc;

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

%% Get edge lengths, angles, and curvatures
lengths     = getEdgeLengths(vert_u(:, 1:3, end_incr), edges);
lengths_p   = getEdgeLengths(vert_f(:, :, end_incr), edges);
error       = (lengths - lengths_p)./lengths*100; % this is backwards. tension is (-)

figure();
plot(1:length(edges), error, 'o--')
xlabel('Edge index')
ylabel('Length error (percent)')

angles_f = foldedCreaseAngles_fast(vert_f(:, :, end_incr), vert_u(:, :, end_incr), edges, adj_faces);
angles_u = foldedCreaseAngles_fast(vert_u(:, :, end_incr), vert_f(:, :, end_incr), edges, adj_faces);
angles = angles_f - angles_u + pi;

curv = getCurv2(vert_u(:, :, end_incr), vert_f(:, :, end_incr), edges, adj_faces, faces);

surf_strain = error./100 + curv.*t/2;
surf_strain(isnan(surf_strain)) = [];
max_strain = max(abs(surf_strain));
mean_strain = mean(abs(surf_strain));
std_strain = std(abs(surf_strain));

%% plotting

if plot_on
    % plot deployed
    deployed = figure('Color', [1 1 1]);
    hold on
    for i = 0:(N-1)
        deployed = plot3dNodesEdges((rot^i*vert_u(:, :, end_incr)'), edges, angles, deployed);
        patch('faces',faces(:,1:3),'vertices',(rot^i*vert_u(:, :, end_incr)')', ...
            'facecolor',[0.7 0.7 0.7], 'facealpha', 0.5, ...
            'edgecolor',[0.3,0.3,0.3], 'edgealpha', 0) ;
    end
    hold off
    axis equal; axis tight; axis off
    
    ax = gca; ax.Clipping = 'off';
    
    % plot folded
    stowed = figure('Color', [1 1 1]);
    inner = [];
    for i = 0:(N-1)
        stowed = plot3dNodesEdges((rot^i*vert_f(:, 1:3, end_incr)'), edges, angles, stowed);
        hold on
        patch('faces',faces(:,1:3),'vertices',(rot^i*vert_f(:, :, end_incr)')', ...
            'facecolor',[0.7 0.7 0.7], 'facealpha', 0.5, ...
            'edgecolor',[0.3,0.3,0.3], 'edgealpha', 0) ;
    end
    hold off
    axis equal; axis tight; axis off
    
    ax = gca; ax.Clipping = 'off';
end

%% Save major fold lines

if save_on
    major_v_u = flip(vert_u([labels.valley_idx labels.brim_valley_idx], :, end_incr), 1);
    major_v_f = flip(vert_f([labels.valley_idx labels.brim_valley_idx], :, end_incr), 1);

    major_m_u = vert_u([labels.mount_idx, labels.brim_mountain_idx], :, end_incr);
    major_m_f = vert_f([labels.mount_idx, labels.brim_mountain_idx], :, end_incr);
    
    save(fullfile(save_path, sprintf("major_folds_c%d_n%d_N%d_rib%d_gamma%d_R%d.mat", [round(c*1000), n, N, rib_d*1000, geo.gamma*100, round(R*100)])), 'major_v_u', 'major_v_f', 'major_m_u', 'major_m_f');
    writematrix(major_v_u(:, 1:2), fullfile(save_path, sprintf('major_v_u_c%d_n%d_N%d_rib%d_gamma%d_R%d.csv', [round(c*1000), n, N, rib_d*1000, geo.gamma*100, round(R*100)])));
    writematrix(major_m_u(:, 1:2), fullfile(save_path, sprintf('major_m_u_c%d_n%d_N%d_rib%d_gamma%d_R%d.csv', [round(c*1000), n, N, rib_d*1000, geo.gamma*100, round(R*100)])));
    
    inner_idxs = vert_u(labels.stat_idxs([1, end:-1:2]), :, end_incr);
    hexagon = inner_idxs;
    
    for i = 1:(N-1)
        hexagon = [hexagon; (rot^i*inner_idxs')'];
    end
    
    hexagon = [hexagon; hexagon(1, :)];
    
    writematrix(hexagon(:, 1:2), fullfile(save_path, sprintf('inner_hexagon_c%d_n%d_N%d_rib%d_gamma%d_R%d.csv', [round(c*1000), n, N, rib_d*1000, geo.gamma*100, round(R*100)])));
    %save(fullfile(save_path, sprintf("101525_converge_c%d_n%d_N%d_rib%d_gamma%d.mat", [round(c*1000), n, N, rib_d*1000, geo.gamma*100])));
end

if save_stability
    nodes_f = vert_f(:, :, end_incr);
    nodes_u = vert_u(:, :, end_incr);
    save(fullfile(save_path, sprintf("040126_stability_c%d_n%d_N%d_rib%d_gamma%d_R%d.mat", [round(c*1000), n, N, rib_d*1000, geo.gamma*100, round(R*100)])), 'nodes_f', 'nodes_u', 'surf_strain', 'error', 'curv', 'E_ax', 'E_cr', 'E_v', 'edges', 'adj_faces', 'faces', 'runtime', 'geo', 'labels');
end

%% Generate unified mesh (This takes a while)
% [unfolded, edges_one, faces_one] = makeFullMesh(vert_u(:, :, end_incr), edges, faces, rot, N);
% [folded, ~, ~] = makeFullMesh(vert_f(:, :, end_incr), edges, faces, rot, N);
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

%% functions
function curv = getCurv2(vert_u, vert_f, edges, adj_faces, faces)
    angles_f = foldedCreaseAngles_fast(vert_f, vert_u, edges, adj_faces);
    angles_u = foldedCreaseAngles_fast(vert_u, vert_f, edges, adj_faces);
    angles_f(abs(angles_f-pi)>pi/4) = nan; % ignore edges connecting ribs

    nEdges = size(edges, 1);    
    h_f = nan(nEdges, 1);
    h_u = nan(nEdges, 1);
    
    for i = 1:nEdges
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);
        
        faces_with_node1        = (faces(:, 1) == p1_index) | (faces(:, 2) == p1_index) | (faces(:, 3) == p1_index);
        faces_with_node2        = (faces(:, 1) == p2_index) | (faces(:, 2) == p2_index) | (faces(:, 3) == p2_index);
        faces_adjacent_to_edge  = faces((faces_with_node1 & faces_with_node2), :);
        
        if size(faces_adjacent_to_edge, 1) == 2
            face1 = faces_adjacent_to_edge(1, :);
            face2 = faces_adjacent_to_edge(2, :);
            
            p3_1_index = face1((face1 ~= p1_index) & (face1 ~= p2_index));
            p3_2_index = face2((face2 ~= p1_index) & (face2 ~= p2_index));

            p1 = vert_f(p1_index, 1:3);
            p2 = vert_f(p2_index, 1:3);
            p1_u = vert_u(p1_index, 1:3);
            p2_u = vert_u(p2_index, 1:3);
                        
            p3_1 = vert_f(p3_1_index, 1:3);
            p3_2 = vert_f(p3_2_index, 1:3); 
            p3_1_u = vert_u(p3_1_index, 1:3);
            p3_2_u = vert_u(p3_2_index, 1:3);

            h1_f = norm(cross(p1-p3_1, p2-p3_1))/norm(p2-p1);
            h2_f = norm(cross(p1-p3_2, p2-p3_2))/norm(p2-p1);

            h1_u = norm(cross(p1_u-p3_1_u, p2_u-p3_1_u))/norm(p2_u-p1_u);
            h2_u = norm(cross(p1_u-p3_2_u, p2_u-p3_2_u))/norm(p2_u-p1_u);

            h_f(i) = (h1_f + h2_f)./2;
            h_u(i) = (h1_u + h2_u)./2;
        end
    end

    curv = (angles_f-pi)./h_f - (angles_u-pi)./h_u;
end

