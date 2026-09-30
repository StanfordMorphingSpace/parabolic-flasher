clear; clc;

set(groot,'defaulttextinterpreter','tex'); 
set(groot,'defaultAxesTickLabelInterpreter','tex'); 
set(groot,'defaultLegendInterpreter','tex');

%addpath('D:\Curved_crease_antennas\MananFunctions')
addpath(genpath('D:\Curved_crease_antennas\MESH2D'))

% add all subfolders
addpath(genpath(pwd))

tip = 0;
%%

%==========================================================================
% Read Data
%==========================================================================
% load matlab data
load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c2117_n50_N8_rib0_gamma10.mat"); % dish (a)
%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\new\040126_stability_c2117_n50_N8_rib0_gamma80.mat"); % dish (a)
%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c350_n50_N10_rib0_gamma20.mat");
%load("D:\Curved_crease_antennas\Reflector\designs\101525_converge_c250_n50_N8_rib0.mat"); 
%load("D:\Curved_crease_antennas\journal_2026\fold\101525_converge_c2117_n50_N8_rib0.mat"); % dish (a) for abaqus comparison

%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c0_n60_N6_rib0_gamma70.mat"); % dish (b)
%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c2500_n80_N16_rib0_gamma70.mat"); % dish (d)
%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c250_n100_N10_rib0_gamma80.mat"); % dish (c)

%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c2117_n50_N8_rib0_gamma10.mat");
%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c2117_n100_N8_rib0_gamma80.mat");
%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\new\040126_stability_c2117_n50_N8_rib0_gamma80.mat");

%load("D:\Curved_crease_antennas\SciTech_2027\fold_pattern\040126_stability_c2117_n50_N6_rib0_gamma10_R25.mat"); % N6 new polycarb

%load("D:\Curved_crease_antennas\SciTech_2027\fold_pattern\040126_stability_c250_n45_N8_rib0_gamma10_R22.mat"); % lunar reflector

%load("D:\Curved_crease_antennas\journal_2026\convergence_comp\040126_stability_c2117_n100_N8_rib0_gamma80");

%brims
%load("D:\Curved_crease_antennas\SciTech_2027\fold_pattern\brims\fold\040126_stability_c2117_n60_N8_rib0_gamma10_R25.mat"); 
%load("D:\Curved_crease_antennas\SciTech_2027\fold_pattern\brims\no_fold\040126_stability_c2117_n60_N8_rib0_gamma10_R25.mat");

vert_u = nodes_u;
vert_f = nodes_f;
%error = -error./100;

% vert_u = vert_u(:, :, end_incr);
% vert_f = vert_f(:, :, end_incr);
%%
lengths     = getEdgeLengths(vert_u, edges);
lengths_p   = getEdgeLengths(vert_f, edges);
error       = (lengths_p - lengths)./lengths;

min_val = min(error');
max_val = max(error');

c_error = round(255*((error)-(min_val))./((max_val)-(min_val)))+1;

cmap = buildDivergingColormap_v7(min_val, max_val);


%% matlab figure
nEdges = size(edges,1);
   
figure('Units', 'normalized', 'Position', [0.05 0.01 0.8, 0.8], 'Color', [1 1 1]);
%figure('Color', [1 1 1]);
%subplot(2, 3, 1)
tile = tiledlayout(1, 2);
nexttile
hold on
for i = 1:nEdges        
    color = cmap(c_error(i), :); % + angles(i)/pi()*[0 1 1];
    
    p1 = vert_u(edges(i,1), :);
    p2 = vert_u(edges(i,2), :);
    if tip
        lh = plot(-[p1(1) p2(1)],-[p1(2) p2(2)], 'Color', color, 'LineWidth', 3.5);
    else
        lh = plot(-[p1(1) p2(1)],-[p1(2) p2(2)], 'Color', color, 'LineWidth', 1.25);
    end
end
   
hold off

axis equal
axis tight
axis off
% view(270, 90)

ax = gca;
%ax.Clipping = 'off';
colormap(ax, cmap);
set(gca,'CLim',[min_val max_val]);
%set(gca,'ColorScale','log')
cb = colorbar;
xlabel(ax, 'Abaqus Min Principal Strain', 'FontSize', 16)
ax.XLabel.Visible = 'on';
ylabel(cb,'Midplane Strain','FontSize',16,'Rotation',270)
cb.Ruler.Exponent = -2;
if tip
    ylim([-0.06 0.001])
end

%% Curvature matlab processing

%curv = -getCurv2(vert_u, vert_f, edges, adj_faces, faces);

min_val = min(curv');
max_val = max(curv');

c_curv = round(255*(curv'-min_val)./(max_val-min_val))+1;
% c_val1 = round(255*(val1-min_val)./(max_val-min_val))+1;
% c_val2 = round(255*(val2-min_val)./(max_val-min_val))+1;

cmap = buildDivergingColormap_v7(min_val, max_val);

%% matlab figure
nEdges = size(edges,1);
   
%subplot(2, 3, 4)
nexttile
hold on
for i = 1:nEdges   
    if isnan(c_curv(i))
        color = [0,0,0];
    else
        color = cmap(c_curv(i), :); % + angles(i)/pi()*[0 1 1];
    end
    
    p1 = vert_u(edges(i,1), :);
    p2 = vert_u(edges(i,2), :);
    if tip
        lh = plot(-[p1(1) p2(1)],-[p1(2) p2(2)], 'Color', color, 'LineWidth', 3.5);
    else
        lh = plot(-[p1(1) p2(1)],-[p1(2) p2(2)], 'Color', color, 'LineWidth', 1.25);
    end
end
   
hold off

axis equal
axis tight
axis off
% view(270, 90)

ax = gca;
%ax.Clipping = 'off';
colormap(ax, cmap);
set(gca,'CLim',[min_val max_val]);
%set(gca,'ColorScale','log')
cb = colorbar;
xlabel(ax, 'Abaqus Min Principal Strain', 'FontSize', 16)
ax.XLabel.Visible = 'on';
ylabel(cb,'Curvature Change, m^{-1}','FontSize',16,'Rotation',270)
if tip
    ylim([-0.06 0.001])
end

%% Surface strain
%surf_strain = max(abs(error + curv.*t/2), abs(error - curv.*t/2));
curv(isnan(curv)) = 0;
surf_strain = [error + curv.*geo.t/2; error - curv.*geo.t/2];

max_strain = max(abs(surf_strain));
mean_strain = mean(abs(surf_strain));
std_strain = std(abs(surf_strain));

ordered_strain = sort(abs(surf_strain));
max_strain_98 = ordered_strain(round(0.98*length(ordered_strain)));

% figure;
% histogram(surf_strain*100, 50);
% xlabel('Surface Strain (%)')
% 
% figure;
% histogram(error*100, 100);
% xlabel('Midplane Strain (%)')

fig = figure('Units', 'inches', 'Position', [1, 1, 8.5, 4], 'Color', 'w');
tlo = tiledlayout(2, 2, 'Padding', 'tight', 'TileSpacing', 'compact');
nexttile(tlo)
% [f, xi] = ksdensity(error*100);
% % Calculate bin width (assumes constant bin width)
% % binWidth = (max(error*100) - min(error*100)) / 50; % Using 50 bins
% % % Scale the curve
% % f_prob = f * binWidth * (length(error*100)/sum(error*100)); 
% plot(xi, f)
histogram(error*100, 50, 'Normalization', 'probability', 'HandleVisibility','off');
xline(mean(error*100), 'r')
xlabel('\epsilon (%)')
legend(sprintf('mean = %.2e', mean(error*100)), 'Location', 'northeast')
ylim([-0.08 1])
xlim([-2, 2])
ylabel('Probability')
fontname('Calibri')
set(gca, 'FontSize', 11)

nexttile(tlo)
histogram((error + curv.*geo.t/2).*100, 50, 'Normalization', 'probability', 'HandleVisibility','off');
xline(mean((error + curv.*geo.t/2).*100), 'r')
xlabel('\epsilon + \kappa s / 2 (%)')
legend(sprintf('mean = %.2e', mean((error + curv.*geo.t/2).*100)), 'Location', 'northeast')
ylim([-0.01 0.15])
xlim([-2, 2])
fontname('Calibri')
set(gca, 'FontSize', 11)

nexttile(tlo)
histogram(curv.*geo.t/2.*100, 50, 'Normalization', 'probability', 'HandleVisibility','off');
xline(mean((curv.*geo.t/2.*100)), 'r')
xlabel('\kappa s /2 (%)')
legend(sprintf('mean = %.2e', mean(curv.*geo.t/2.*100)), 'Location', 'northeast')
ylim([-0.01 0.15])
xlim([-2, 2])
ylabel('Probability')
fontname('Calibri')
set(gca, 'FontSize', 11)

nexttile(tlo)
histogram((error - curv.*geo.t/2).*100, 50, 'Normalization', 'probability', 'HandleVisibility','off');
xline(mean((error - curv.*geo.t/2).*100), 'r')
xlabel('\epsilon - \kappa s / 2 (%)')
legend(sprintf('mean = %.2e', mean((error - curv.*geo.t/2).*100)), 'Location', 'northeast')
ylim([-0.01 0.15])
xlim([-2, 2])
fontname('Calibri')
set(gca, 'FontSize', 11)

%% plot energies
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

%% plot shapes

%% Get edge lengths, angles, and curvatures

angles_f = foldedCreaseAngles_fast(vert_f(:, :), vert_u(:, :), edges, adj_faces);
angles_u = foldedCreaseAngles_fast(vert_u(:, :), vert_f(:, :), edges, adj_faces);
angles = angles_f - angles_u + pi;

beta    = 2*pi()/geo.N;
rot     = [ cos(beta), -sin(beta), 0;...
            sin(beta), cos(beta), 0;...
            0, 0, 1];
% 
% % plot deployed
% deployed = figure('Color', [1 1 1]);
% hold on
% for i = 1 %0:(geo.N-1)
%     deployed = plot3dNodesEdges((rot^i*vert_u(:, :)'), edges, angles, deployed);
%     %patch('faces',faces(:,1:3),'vertices',(rot^i*vert_u(:, :)')', ...
%     %    'facecolor',[0.7 0.7 0.7], 'facealpha', 0.5, ...
%     %    'edgecolor',[0.3,0.3,0.3], 'edgealpha', 0) ;
% end
% hold off
% axis equal; axis tight; axis off
% 
% ax = gca; ax.Clipping = 'off';
% % 
% % % plot folded
% % stowed = figure('Color', [1 1 1]);
% % inner = [];
% % for i = 0:(geo.N-1)
% %     stowed = plot3dNodesEdges((rot^i*vert_f(:, 1:3)'), edges, angles, stowed);
% %     hold on
% %     patch('faces',faces(:,1:3),'vertices',(rot^i*vert_f(:, :)')', ...
% %         'facecolor',[0.7 0.7 0.7], 'facealpha', 0.5, ...
% %         'edgecolor',[0.3,0.3,0.3], 'edgealpha', 0) ;
% % end
% % hold off
% % axis equal; axis tight; axis off
% % 
% % ax = gca; ax.Clipping = 'off';

%% Generate unified mesh (This takes a while)
beta    = 2*pi()/geo.N;

rot     = [ cos(beta), -sin(beta), 0;...
            sin(beta), cos(beta), 0;...
            0, 0, 1];
angles = foldedCreaseAngles_fast(vert_f, vert_u, edges, adj_faces);
[unfolded, folded, edges_one, faces_one, angles_one] = makeFullMesh_v2(vert_u, vert_f, edges, faces, angles, rot, geo.N);

%%
folded_one = figure('Color', [1 1 1]);
patch('faces',faces_one,'vertices',folded, ...
        'facecolor',[0.7 0.7 0.7], 'facealpha', 0.4, ...
        'edgecolor',[1,0,0], 'edgealpha', 0) ;

folded_one = plot3dNodesEdges(folded', edges_one, angles_one, folded_one);

axis equal; axis tight; axis off; view(-37.5, 30)

ax = gca; ax.Clipping = 'off';

deployed_one = figure('Color', [1 1 1]);
% patch('faces',faces_one,'vertices',unfolded, ...
%         'facecolor',[0.7 0.7 0.7], 'facealpha', 0.5, ...
%         'edgecolor',[1,0,0], 'edgealpha', 0) ;

deployed_one = plot3dNodesEdges(unfolded', edges_one, angles_one, deployed_one);

axis equal; axis tight; axis off; view(-37.5, 30)

ax = gca; ax.Clipping = 'off';

%%

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

function R = alignZAxisWithVector(d)
    % Ensure the input vector is a column vector
    d = d(:);
    
    % Step 1: Normalize the desired direction vector d to create the new z-axis
    z_prime = d / norm(d);
    
    % Step 2: Choose a vector a that is not parallel to z_prime
    % Here, we use [1; 0; 0] as long as it's not parallel to z_prime
    if abs(z_prime(1)) < 1 - 1e-10  % Check if z_prime is not [±1; 0; 0]
        a = [1; 0; 0];
    else
        a = [0; 1; 0]; % Use [0; 1; 0] if z_prime is [±1; 0; 0]
    end
    
    % Step 3: Compute the new y-axis by taking the cross product
    y_prime = cross(z_prime, a);
    y_prime = y_prime / norm(y_prime);  % Normalize the new y-axis
    
    % Step 4: Compute the new x-axis as the cross product of y_prime and z_prime
    x_prime = cross(y_prime, z_prime);
    
    % Step 5: Form the transformation matrix R
    R = [x_prime, y_prime, z_prime];
end

function cmap = buildDivergingColormap_v7(minVal, maxVal)
    % Diverging colormap: deep red -> orange -> light orange -> white ->
    % light blue -> vibrant blue -> dark purple, with the lightest inner
    % colors always exactly 6 indices away from white (zero).
    colors = [ ...
        0.75, 0.06, 0.06;   % Deep red (Far left)
        0.76, 0.35, 0.0;    % orange (Intermediate 1 left)
        0.9,  0.74, 0.60;   % light orange (Intermediate 2 left)
        1.00, 1.00, 1.00;   % Pure White (Center/Zero)
        0.60, 0.70, 0.90;   % Light Sky Blue (Intermediate 1 right)
        0.05, 0.25, 0.60;   % Vibrant Blue (Intermediate 2 right)
        0.39, 0.0,  0.76];  % Dark Purple (Far right)

    if minVal == maxVal
        cmap = repmat(colors(4, :), 256, 1);
        return;
    end

    f_zero = (0 - minVal) / (maxVal - minVal);

    if f_zero > 0 && f_zero < 1
        % --- ZERO IS IN RANGE ---
        zero_idx = round(1 + 255 * f_zero);

        idx0 = zeros(1, 7);
        idx0(1) = 1;
        idx0(7) = 256;
        idx0(4) = zero_idx;
        idx0(3) = zero_idx - 8;
        idx0(5) = zero_idx + 8;
        idx0(2) = round(1 + 0.7 * (idx0(3) - 1));
        idx0(6) = round(idx0(5) + 0.3 * (256 - idx0(5)));

        colors_to_use = colors;
        idx0 = max(1, min(256, idx0));

        keep_mask = true(1, 7);
        for i = 3:-1:1
            if idx0(i) >= idx0(i+1)
                keep_mask(i) = false;
                idx0(i) = idx0(i+1);
            end
        end
        for i = 5:7
            if idx0(i) <= idx0(i-1)
                keep_mask(i) = false;
                idx0(i) = idx0(i-1);
            end
        end

    elseif f_zero <= 0
        % --- ZERO IS NOT IN RANGE (data entirely positive) ---
        idx0 = [1, 1 + 6, round((1 + 6) + 0.5 * (256 - (1 + 6))), 256];
        colors_to_use = colors(4:7, :);

        keep_mask = true(1, 4);
        for i = 2:4
            if idx0(i) <= idx0(i-1)
                keep_mask(i) = false;
                idx0(i) = idx0(i-1);
            end
        end

    else
        % --- ZERO IS NOT IN RANGE (data entirely negative) ---
        idx0 = [1, round(1 + 0.5 * ((256 - 6) - 1)), 256 - 6, 256];
        colors_to_use = colors(1:4, :);

        keep_mask = true(1, 4);
        for i = 3:-1:1
            if idx0(i) >= idx0(i+1)
                keep_mask(i) = false;
                idx0(i) = idx0(i+1);
            end
        end
    end

    idx0_unique = idx0(keep_mask);
    colors_unique = colors_to_use(keep_mask, :);

    x = (idx0_unique - 1) / 255;
    xq = linspace(0, 1, 256);

    if length(idx0_unique) >= 3
        method = 'pchip';
    else
        method = 'linear';
    end

    r = interp1(x, colors_unique(:,1), xq, method);
    g = interp1(x, colors_unique(:,2), xq, method);
    b = interp1(x, colors_unique(:,3), xq, method);

    cmap = max(0, min(1, [r' g' b']));
end

