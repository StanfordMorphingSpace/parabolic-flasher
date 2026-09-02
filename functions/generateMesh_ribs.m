function [vert_u, vert_f, vert_ref, labels, edges, faces] = generateMesh_ribs(A, N, h, n, R, surf_func, plot_on, rib_d, brim_R, brim_func)
% Builds the unfolded (projected onto z = surf_func(r), r = sqrt(x^2+y^2)),
% folded, and flat reference meshes for one valley/mountain gore pair of an
% N-gore flasher, and optionally appends stiffening ribs along the valley
% fold and/or a brim (projected onto z = brim_func(r)) beyond radius R.
%
% labels is a struct of node-index groups used by the Jacobian/b constraint
% functions: stat_idxs, outer_idx, cone_idx, valley_idx, valley_rot_idx
% (the rotated copy of valley_idx, same order, for the sector-periodicity
% constraint), mount_idx, and (only when rib_d ~= 0) ribs_idx and
% ribs_lower_idx (the bottom row of rib nodes). When brim_R ~= 0, also:
% brim_idx (all newly added brim nodes), brim_outer_idx (brim nodes at
% radius brim_R), brim_valley_idx/brim_valley_rot_idx/brim_mountain_idx
% (continuations of valley_idx/valley_rot_idx/mount_idx beyond R, same
% ordering convention, excluding the r=R seam already covered by those).

%% Sector geometry setup
beta    = 2*pi()/N;
delta_z = 2*A*sin(beta/2)*tan(beta/2); % change in height of mountain fold for a zero-thickness flasher

rot     = [ cos(beta), -sin(beta), 0;...
            sin(beta), cos(beta), 0;...
            0, 0, 1];
rot_gen  = @(theta) [ cos(theta), -sin(theta), 0;...
            sin(theta), cos(theta), 0;...
            0, 0, 1];

% endpoint where a line through p1,p2 crosses a circle of the given radius
line_circle_x = @(p1, p2, rad) ((p1(2)-p2(2))*(p2(1)*p1(2) - p1(1)*p2(2)) + sqrt((p1(1)-p2(1))^2*(rad^2*norm(p1-p2)^2 - (p2(1)*p1(2) - p1(1)*p2(2))^2)))/(norm(p1-p2)^2);
line_circle_y = @(p1, p2, rad) (-p2(1)^3*p1(2) + p1(1)^3*p2(2) + p1(1)*p2(1)^2*(2*p1(2) + p2(2)) - p1(1)^2*p2(1)*(p1(2) + 2*p2(2)) + (p1(2)-p2(2))*sqrt((p1(1)-p2(1))^2*(rad^2*norm(p1-p2)^2 - (p2(1)*p1(2) - p1(1)*p2(2))^2)))/((p1(1)-p2(1))*norm(p1-p2)^2);

%% Valley (v) and mountain (m) fold line endpoints across one gore sector
init_vec = [linspace(A, R, 2); zeros(1, 2); zeros(1, 2)];
v_init = rot_gen((pi-beta)/2)*(init_vec-init_vec(:, 1))+[A;0;0];
m_init = rot_gen(pi/2)*(init_vec-init_vec(:, 1))+[A;0;0];

v_end = [line_circle_x(v_init(:, 1), v_init(:, 2), R); line_circle_y(v_init(:, 1), v_init(:, 2), R); 0];
m_end = [m_init(1, 1), sqrt(R^2-m_init(1, 1)^2)];

v_plane = [linspace(A, v_end(1), n); linspace(0, v_end(2), n); zeros(1, n)];
m_plane = [linspace(A, m_end(1), n); linspace(0, m_end(2), n); zeros(1, n)];

% find the last subdivision near the tip closer together than one mesh spacing,
% so it can be triangulated by hand instead of by refine2
idx = 0;
dis = 0;
while dis <= (R-A)/n
    idx = idx+1;
    dis = norm(v_plane(:, idx)-m_plane(:, idx));
end
idx = idx-1;

hfun = (R-A)/n; % target triangle edge length for refine2

%% Gore 1: region between the valley and mountain fold lines
end_edge_d = norm(v_plane(:, end)-m_plane(:, end));
end_n = round(end_edge_d/((R-A)/n))+1;
end_theta = angleBetweenVectors(v_plane(:, end)', m_plane(:, end)');
end_thetas = linspace(0, end_theta, end_n);
end_edge = [];
for i = 2:(end_n-1)
    end_edge = [end_edge, rot_gen(-end_thetas(i))*m_plane(:, end)];
end
end_edge_gore1 = end_edge; % gore 2 reuses the name "end_edge" for its own (different) arc below

nodes1 = [flip(v_plane(:, idx:end), 2), m_plane(:, idx:end), end_edge];
gore1 = makeGoreLoop(length(nodes1));

[vert1, tria1] = triangulateGore(nodes1, gore1, hfun);

% add back the tip triangles refine2 trims for being smaller than hfun
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

%% Gore 2: mirrored region on the other side of the mountain fold
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
gore2 = makeGoreLoop(length(nodes2));

[vert2, tria2] = triangulateGore(nodes2, gore2, hfun);

% add back the tip triangles refine2 trims for being smaller than hfun
if idx > 2
    num_vert = size(vert2, 1);
    edge_m = [flip(inner_edge, 2), rot*v_plane];
    vert2 = [vert2;  m_plane(1:2, 1:(idx-1))'; edge_m(1:2, 1:(idx-2))'];
    tip_idx2 = find(ismember(vert2, edge_m(1:2, idx-1)', 'rows'));
    tip_idx1 = find(ismember(vert2, m_plane(1:2, idx)', 'rows'));

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

    tria2 = [tria2; num_vert+1, tip_idx1, tip_idx2];
end

%% Assign mesh vertices to ruling triangles
% Ruling triangles fan out from the fold lines to the tip/outer edge; every
% mesh vertex is either a ruling-line node itself, or lies inside exactly
% one ruling triangle and gets barycentrically interpolated within it later.
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

vert1 = [vert1, zeros(length(vert1), 1)];
labels1 = assignRulingLabels(vert1, ruling_nodes1, ruling_tris1);

vert2 = [vert2, zeros(length(vert2), 1)];
labels2 = assignRulingLabels(vert2, ruling_nodes2, ruling_tris2);

%% Project unfolded mesh onto the surface z = surf_func(r), r = sqrt(x^2+y^2)
vert1_p = vert1;
vert1_p(:, 3) = surf_func(sqrt(vert1_p(:, 1).^2 + vert1_p(:, 2).^2));

vert2_p = vert2;
vert2_p(:, 3) = surf_func(sqrt(vert2_p(:, 1).^2 + vert2_p(:, 2).^2));

v_plane_p = v_plane;
v_plane_p(3, :) = surf_func(sqrt(v_plane_p(1, :).^2 + v_plane_p(2, :).^2));

m_plane_p = m_plane;
m_plane_p(3, :) = surf_func(sqrt(m_plane_p(1, :).^2 + m_plane_p(2, :).^2));

m_plane_ext_p = m_plane_ext;
m_plane_ext_p(3, :) = surf_func(sqrt(m_plane_ext_p(1, :).^2 + m_plane_ext_p(2, :).^2));

inner_edge_p = inner_edge;
inner_edge_p(3, :) = surf_func(sqrt(inner_edge_p(1, :).^2 + inner_edge_p(2, :).^2));

ruling_nodes1_p = [v_plane_p, m_plane_p(:, 2:end)];
ruling_nodes2_p = [m_plane_ext_p, flip(inner_edge_p, 2), rot*v_plane_p(:, 1:end)];

m_folded    = nan(3, length(v_plane_p));
v_folded    = nan(3, length(v_plane_p));
v_folded(:, 1) = [A; 0; min(v_plane_p(3, :))];
m_folded(:, 1) = v_folded(:, 1);

%% Folded initial condition
% Walk each ruling line out along a spiral of constant edge length (matching
% the projected mesh's edge lengths) to get a valid folded starting shape.
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

ruling_nodes1_f = [v_folded, m_folded(:, 2:end)];
ruling_nodes2_f = [m_folded_ext, flip(inner_edge_p, 2), rot*v_folded(:, 1:end)];

cone_idx = sum(theta_m < beta); % nodes within one sector of the tip; not constrained to the spiral

%% Place subdivided mesh vertices within their ruling triangles
% Vertices on a ruling line take that line's folded position directly;
% the rest are barycentrically interpolated within their enclosing triangle.
vert1_f = interpolateFoldedVertices(vert1_p, labels1, ruling_nodes1_p, ruling_tris1, ruling_nodes1_f);
vert2_f = interpolateFoldedVertices(vert2_p, labels2, ruling_nodes2_p, ruling_tris2, ruling_nodes2_f);

%% Combine both gores into one section, merging shared vertices along the mountain fold
vert_ref = vert1;
vert_u = vert1_p;
vert_f = vert1_f;
tria_temp = zeros(size(tria2));
% raw_labels encodes, per node: k in [1,n] on the valley line, n+k on the
% mountain line (k>1), -k on the rotated copy of valley node k seen from the
% far side of gore 2, or nan elsewhere. Collapsed into named index groups below.
raw_labels = labels1(:, 2);

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
            raw_labels = [raw_labels; -ruling_idx];
        else
            raw_labels = [raw_labels; nan];
        end
    else
        tria_temp(tria_idx) = found_idx;
    end

end

faces = [tria1; tria_temp];

% derive the (undirected, deduplicated) edge list from the assembled faces
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

%% Find indices of constrained (fixed) nodes
% The tip and the inner edge between the two gores stay put during dynamic
% relaxation; vert_f is snapped back to vert_u there to enforce that exactly.
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

%% Collapse raw_labels into named node-index groups
% valley_idx(k) and valley_rot_idx(k) refer to the same ruling point (the
% k-th subdivision along the valley line, and its rotated copy across gore
% 2), kept in matching order so the rotational-periodicity constraint can
% pair them up directly instead of searching for -label(x) each time.
valley_idx = zeros(1, n);
valley_rot_idx = zeros(1, n);
for k = 1:n
    valley_idx(k) = find(raw_labels == k, 1);
    valley_rot_idx(k) = find(raw_labels == -k, 1);
end
% mount_idx(k) is the (k+1)-th subdivision along the mountain line (there is
% no separate node for k=1: it's shared with valley_idx(1)), ordered outward
mount_idx = zeros(1, n-1);
for k = 1:(n-1)
    mount_idx(k) = find(raw_labels == n+k, 1);
end

%% Ribs (if a non-zero rib depth is specified)
% Extrude a stack of rib rows straight down (in -z) from the valley (v)
% ruling line, connecting each row to the mesh above and the row below it.
if rib_d > 0
    nNodes = length(vert_u);
    rib_n = max(round(rib_d/((R-A)/n))+1, 2);
    rib_z = linspace(0, rib_d, rib_n);

    % unfolded
    rib_nodes_u = zeros(3, length(v_plane_p)*(rib_n-1));
    for j = 1:(rib_n-1)
        rib_nodes_u(:, ((j-1)*length(v_plane_p)+1):(j*length(v_plane_p))) = v_plane_p - [zeros(2, length(v_plane_p)); rib_z(j+1).*ones(1, length(v_plane_p))];
    end

    % folded
    rib_nodes_f = zeros(3, length(v_folded)*(rib_n-1));
    for j = 1:(rib_n-1)
        rib_nodes_f(:, ((j-1)*length(v_folded)+1):(j*length(v_folded))) = v_folded - [zeros(2, length(v_folded)); rib_z(j+1).*ones(1, length(v_folded))];
    end

    % flat reference
    rib_nodes_ref = zeros(3, length(v_plane)*(rib_n-1));
    for j = 1:(rib_n-1)
        rib_nodes_ref(:, ((j-1)*length(v_plane)+1):(j*length(v_plane))) = v_plane - [zeros(2, length(v_plane)); rib_z(j+1).*ones(1, length(v_plane))];
    end

    rib_edges = [];
    rib_faces = [];
    for k = 1:(n-1) % first row, connecting the mesh to the top of the ribs
        idx1 = valley_idx(k);
        idx2 = valley_idx(k+1);
        rib_edges = [rib_edges; idx1, nNodes+k; idx1, nNodes+k+1];
        rib_faces = [rib_faces; idx1, idx2, nNodes+k+1; idx1, nNodes+k+1, nNodes+k];
    end
    rib_edges = [rib_edges; valley_idx(n), nNodes+n];
    for j = 1:(rib_n-2) % remaining rows, connecting each rib row to the next
        for k = 1:(n-1)
            base_idx = nNodes + (j-1)*n;
            rib_edges = [rib_edges; base_idx+k, base_idx+k+1; base_idx+k, base_idx+k+n; base_idx+k, base_idx+k+n+1];
            rib_faces = [rib_faces; base_idx+k, base_idx+k+1, base_idx+k+n+1; base_idx+k, base_idx+k+n+1, base_idx+k+n];
        end
        rib_edges = [rib_edges; nNodes+j*n, nNodes+(j+1)*n]; % end edge
    end
    for k = 1:(n-1)
        rib_edges = [rib_edges; nNodes+(rib_n-2)*n+k, nNodes+(rib_n-2)*n+k+1]; % bottom edge
    end

    vert_u = [vert_u; rib_nodes_u'];
    vert_f = [vert_f; rib_nodes_f'];
    vert_ref = [vert_ref; rib_nodes_ref'];
    edges = [edges; rib_edges];
    faces = [faces; rib_faces];
    stat_idxs = [stat_idxs, nNodes+1+n.*(0:(rib_n-2))];
    outer_idx = [outer_idx, nNodes+n+n.*(0:(rib_n-2))];

    ribs_idx = (nNodes+1):(nNodes+(rib_n-1)*n); % all newly added rib nodes
    ribs_lower_idx = (nNodes+(rib_n-2)*n+1):(nNodes+(rib_n-1)*n); % the bottom row only
end

%% Brim (if a non-zero brim outer radius is specified)
% A second valley/mountain wedge pair, structured exactly like gore 1/gore 2
% but spanning r = R to r = brim_R instead of A to R (so no tip to snip),
% continuing the same valley/mountain rays and projecting onto brim_func.
% The r=R boundary is shared with the base mesh (not re-triangulated): the
% seam nodes are given the exact existing coordinates so the merge below
% recognizes and reuses them instead of duplicating them.
if brim_R > 0
    nNodes_before_brim = length(vert_u);
    brim_n = max(round((brim_R-R)/((R-A)/n))+1, 2);

    v_brim_end = [line_circle_x(v_init(:, 1), v_init(:, 2), brim_R); line_circle_y(v_init(:, 1), v_init(:, 2), brim_R); 0];
    m_brim_end = [A; sqrt(brim_R^2-A^2); 0];

    v_plane_brim = [linspace(v_end(1), v_brim_end(1), brim_n); linspace(v_end(2), v_brim_end(2), brim_n); zeros(1, brim_n)];
    m_plane_brim = [linspace(m_end(1), m_brim_end(1), brim_n); linspace(m_end(2), m_brim_end(2), brim_n); zeros(1, brim_n)];

    % brim gore 1: valley-side wedge, outer arc sweeps valley->mountain (as gore 2's does)
    end_theta_brim1 = angleBetweenVectors(v_brim_end', m_brim_end');
    end_n_brim1 = round(norm(v_brim_end-m_brim_end)/((R-A)/n))+1;
    end_thetas_brim1 = linspace(0, end_theta_brim1, end_n_brim1);
    outer_arc_brim1 = [];
    for i = 2:(end_n_brim1-1)
        outer_arc_brim1 = [outer_arc_brim1, rot_gen(end_thetas_brim1(i))*v_brim_end];
    end

    nodes1_brim = [v_plane_brim, outer_arc_brim1, flip(m_plane_brim, 2), end_edge_gore1];
    gore1_brim = makeGoreLoop(length(nodes1_brim));
    [vert1_brim, tria1_brim] = triangulateGore(nodes1_brim, gore1_brim, hfun);
    vert1_brim = [vert1_brim, zeros(length(vert1_brim), 1)];

    % brim gore 2: mountain-side wedge, outer arc sweeps mountain->valley (as gore 1's does)
    v_brim_end_rot = rot*v_brim_end;
    end_theta_brim2 = angleBetweenVectors(v_brim_end_rot', m_brim_end');
    end_n_brim2 = round(norm(v_brim_end_rot-m_brim_end)/((R-A)/n))+1;
    end_thetas_brim2 = linspace(0, end_theta_brim2, end_n_brim2);
    outer_arc_brim2 = [];
    for i = 2:(end_n_brim2-1)
        outer_arc_brim2 = [outer_arc_brim2, rot_gen(-end_thetas_brim2(i))*m_brim_end];
    end

    nodes2_brim = [m_plane_brim, outer_arc_brim2, flip(rot*v_plane_brim, 2), end_edge];
    gore2_brim = makeGoreLoop(length(nodes2_brim));
    [vert2_brim, tria2_brim] = triangulateGore(nodes2_brim, gore2_brim, hfun);
    vert2_brim = [vert2_brim, zeros(length(vert2_brim), 1)];

    % brim ruling "ladder": no shared apex (unlike the gore fans), just a
    % strip of quads between the valley/mountain (or mountain/rotated-valley) lines
    ruling_nodes1_brim = [v_plane_brim, m_plane_brim];
    ruling_tris1_brim = makeLadderTris(brim_n);
    ruling_nodes2_brim = [m_plane_brim, rot*v_plane_brim];
    ruling_tris2_brim = makeLadderTris(brim_n);

    labels1_brim = assignRulingLabels(vert1_brim, ruling_nodes1_brim, ruling_tris1_brim);
    labels2_brim = assignRulingLabels(vert2_brim, ruling_nodes2_brim, ruling_tris2_brim);

    % project the brim onto brim_func, but keep every r=R seam node (not just
    % the two corners - the whole reused end_edge/end_edge_gore1 arcs are at
    % r=R too) exactly on the existing surf_func-projected value, since
    % brim_func(R) is not guaranteed bit-identical to surf_func(R)
    v_plane_brim_p = projectBrimNodes(v_plane_brim', R, surf_func, brim_func)';
    m_plane_brim_p = projectBrimNodes(m_plane_brim', R, surf_func, brim_func)';
    vert1_brim_p = projectBrimNodes(vert1_brim, R, surf_func, brim_func);
    vert2_brim_p = projectBrimNodes(vert2_brim, R, surf_func, brim_func);

    ruling_nodes1_brim_p = [v_plane_brim_p, m_plane_brim_p];
    ruling_nodes2_brim_p = [m_plane_brim_p, rot*v_plane_brim_p];

    % continue the folded spiral outward, starting exactly from the existing
    % folded seam position so the folded mesh is a seamless continuation too
    v_folded_brim = zeros(3, brim_n);
    v_folded_brim(:, 1) = v_folded(:, end);
    theta_v_brim = theta_v(end);
    for i = 2:brim_n
        dist_v = norm(v_plane_brim_p(:, i-1)-v_plane_brim_p(:, i));
        fun_v = @(theta) (dist_v - norm(v_folded_brim(:, i-1)-v_spiral(theta)))^2;
        th_v = fmincon(fun_v, theta_v_brim(end)+beta, [], [], [], [], theta_v_brim(end), theta_v_brim(end)+beta, [], options);
        v_folded_brim(:, i) = v_spiral(th_v);
        theta_v_brim = [theta_v_brim, th_v];
    end

    m_folded_brim = zeros(3, brim_n);
    m_folded_brim(:, 1) = m_folded(:, end);
    theta_m_brim = theta_m(end);
    for i = 2:brim_n
        dist_m = norm(m_plane_brim_p(:, i-1)-m_plane_brim_p(:, i));
        fun_m = @(theta) (dist_m - norm(m_folded_brim(:, i-1)-m_spiral(theta)-[0;0;theta*delta_z/beta]))^2;
        th_m = fmincon(fun_m, theta_m_brim(end)+beta, [], [], [], [], theta_m_brim(end), theta_m_brim(end)+beta, [], options);
        m_folded_brim(:, i) = m_spiral(th_m) + [0;0;th_m*delta_z/beta];
        theta_m_brim = [theta_m_brim, th_m];
    end

    ruling_nodes1_brim_f = [v_folded_brim, m_folded_brim];
    ruling_nodes2_brim_f = [m_folded_brim, rot*v_folded_brim];

    vert1_brim_f = interpolateFoldedVertices(vert1_brim_p, labels1_brim, ruling_nodes1_brim_p, ruling_tris1_brim, ruling_nodes1_brim_f);
    vert2_brim_f = interpolateFoldedVertices(vert2_brim_p, labels2_brim, ruling_nodes2_brim_p, ruling_tris2_brim, ruling_nodes2_brim_f);

    % merge into the combined mesh, reusing the shared r=R seam nodes exactly
    [vert_u, vert_f, vert_ref, vert_map1] = mergeMeshPiece(vert_u, vert_f, vert_ref, vert1_brim_p, vert1_brim_f, vert1_brim);
    [vert_u, vert_f, vert_ref, vert_map2] = mergeMeshPiece(vert_u, vert_f, vert_ref, vert2_brim_p, vert2_brim_f, vert2_brim);

    brim_faces = [vert_map1(tria1_brim); vert_map2(tria2_brim)];
    faces = [faces; brim_faces];

    for i = 1:size(brim_faces, 1)
        edge1 = brim_faces(i, 1:2);
        edge2 = brim_faces(i, 2:3);
        edge3 = brim_faces(i, [1 3]);
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

    % named brim groups: exclude ruling index 1 (the r=R seam), already
    % covered by valley_idx(n)/mount_idx(n-1)/valley_rot_idx(n)
    brim_valley_idx = zeros(1, brim_n-1);
    brim_mountain_idx = zeros(1, brim_n-1);
    brim_valley_rot_idx = zeros(1, brim_n-1);
    for k = 2:brim_n
        brim_valley_idx(k-1) = vert_map1(find(labels1_brim(:, 2) == k, 1));
        brim_mountain_idx(k-1) = vert_map1(find(labels1_brim(:, 2) == brim_n+k, 1));
        brim_valley_rot_idx(k-1) = vert_map2(find(labels2_brim(:, 2) == brim_n+k, 1));
    end

    brim_idx = unique([vert_map1(vert_map1 > nNodes_before_brim); vert_map2(vert_map2 > nNodes_before_brim)])';

    brim_outer_idx = [];
    for i = brim_idx
        if ismembertol(norm(vert_u(i, 1:2)), brim_R, 1e-7) && i ~= brim_valley_rot_idx(end)
            brim_outer_idx = [brim_outer_idx, i];
        end
    end
end

%% Collect everything into the labels struct
labels.stat_idxs = stat_idxs;
labels.outer_idx = outer_idx;
labels.cone_idx = cone_idx;
labels.valley_idx = valley_idx;
labels.valley_rot_idx = valley_rot_idx;
labels.mount_idx = mount_idx;
if rib_d > 0
    labels.ribs_idx = ribs_idx;
    labels.ribs_lower_idx = ribs_lower_idx;
end
if brim_R > 0
    labels.brim_idx = brim_idx;
    labels.brim_outer_idx = brim_outer_idx;
    labels.brim_valley_idx = brim_valley_idx;
    labels.brim_valley_rot_idx = brim_valley_rot_idx;
    labels.brim_mountain_idx = brim_mountain_idx;
end

%% Plot unfolded and folded mesh (optional, for visual sanity-checking)
if plot_on
    % unfolded mesh, colored by ruling triangle
    figure('Color', [1 1 1]);
    patch('faces',faces(:,1:3),'vertices',vert_u(:, 1:2), ...
        'facecolor','w', 'facealpha', 0, ...
        'edgecolor',[1,0,0], 'edgealpha', 0.2) ;
    hold on; axis image off;

    patch('faces',ruling_tris1(:,1:2),'vertices', ruling_nodes1(1:2, :)', ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    gscatter(vert1_p(:, 1), vert1_p(:, 2), labels1(:, 1))

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

    if brim_R > 0
        v_verts_brim = v_plane_brim(1:2, :);
        m_verts_brim = m_plane_brim(1:2, :);
        rot_v_verts_brim = (rot*v_plane_brim(:,:));
        plot(v_verts_brim(1, :), v_verts_brim(2, :), 'b--', 'LineWidth', 1.5)
        plot(m_verts_brim(1, :), m_verts_brim(2, :), 'r--', 'LineWidth', 1.5)
        plot(rot_v_verts_brim(1, :), rot_v_verts_brim(2, :), 'g--', 'LineWidth', 1.5)
    end

    inner = [v_plane(1:2, 1), v_plane(1:2, 1)];
    plot(inner(1, :), inner(2, :), 'k-', 'LineWidth', 1.5)
    axis equal
    axis tight
    axis off

    ax = gca;
    ax.Clipping = 'off';

    % folded mesh
    figure('Color', [1 1 1]);
    patch('faces',faces(:,1:3),'vertices',vert_f(:, 1:3), ...
        'facecolor','w', 'facealpha', 0, ...
        'edgecolor',[1,0,0], 'edgealpha', 0.2) ;
    patch('faces',ruling_tris2(:,1:2),'vertices', ruling_nodes2_f(1:3, :)', ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    hold on; axis image off;

    v_verts = v_folded(1:3, :);
    m_verts = m_folded(1:3, :);
    plot3(v_verts(1, :), v_verts(2, :), v_verts(3, :), 'b-', 'LineWidth', 1.5)
    plot3(v_verts(1, :), v_verts(2, :), v_verts(3, :), 'k.', 'MarkerSize', 1)
    plot3(m_verts(1, :), m_verts(2, :), m_verts(3, :), 'r-', 'LineWidth', 1.5)
    plot3(m_verts(1, :), m_verts(2, :), m_verts(3, :), 'k.', 'MarkerSize', 1)

    if brim_R > 0
        v_verts_brim = v_folded_brim(1:3, :);
        m_verts_brim = m_folded_brim(1:3, :);
        rot_v_verts_brim = rot*v_folded_brim(1:3, :);
        plot3(v_verts_brim(1, :), v_verts_brim(2, :), v_verts_brim(3, :), 'b--', 'LineWidth', 1.5)
        plot3(m_verts_brim(1, :), m_verts_brim(2, :), m_verts_brim(3, :), 'r--', 'LineWidth', 1.5)
        plot3(rot_v_verts_brim(1, :), rot_v_verts_brim(2, :), rot_v_verts_brim(3, :), 'g--', 'LineWidth', 1.5)
    end
    axis equal
    axis tight
    axis off

    ax = gca;
    ax.Clipping = 'off';
end

end

function gore = makeGoreLoop(nNodesInGore)
% Closed-polygon boundary edge list (for refine2) over points 1..nNodesInGore in order
gore = [nNodesInGore, 1];
for i = 1:(nNodesInGore-1)
    gore = [gore; i, i+1];
end
end

function [vert, tria] = triangulateGore(nodes, gore, hfun)
% Triangulate a gore's boundary polygon with refine2, then drop the
% degenerate slivers (near-zero interior angle) refine2 sometimes leaves
% along nearly-straight boundary segments.
opts.ref1 = 'preserve';
[vert, ~, tria, ~] = refine2(nodes(1:2, :)', gore, [], opts, hfun);

alphas_test = real(getAlphas([vert, zeros(length(vert), 1)], tria));
bad_idxs = [];
for i = 1:length(alphas_test)
    if abs(alphas_test(i, 1)) < 1e-6 || abs(alphas_test(i, 2)) < 1e-6 || abs(alphas_test(i, 3)) < 1e-6
        bad_idxs = [bad_idxs, i];
    end
end
tria(bad_idxs, :) = [];
end

function labels = assignRulingLabels(vert, ruling_nodes, ruling_tris)
% For each mesh vertex, record which ruling triangle contains it (column 1)
% and, if it coincides exactly with a ruling-line node, that node's index
% (column 2) so it can later be pinned to the folded ruling node directly
% instead of barycentrically interpolated.
labels = zeros(length(vert), 2);
for i = 1:length(vert)
    ruling_idx = find(ismember(ruling_nodes(1:2, :)', vert(i, 1:2), 'rows'));
    if ~isempty(ruling_idx)
        labels(i, 2) = ruling_idx;
    else
        labels(i, 2) = nan;
    end
    for j = 1:length(ruling_tris)
        [in, on] = inpolygon(vert(i, 1), vert(i, 2), ruling_nodes(1, ruling_tris(j, :)), ruling_nodes(2, ruling_tris(j, :)));
        if in == 1
            if on == 0
                labels(i, 1) = j;
            elseif isempty(ruling_idx)
                labels(i, 1) = j;
            else
                labels(i, 1) = nan;
            end
            break
        end
    end
end
labels(labels==0) = max(labels(:, 1)); % bug workaround
end

function vert_f = interpolateFoldedVertices(vert_p, labels, ruling_nodes_p, ruling_tris, ruling_nodes_f)
% Map each projected mesh vertex to its folded position: vertices that
% coincide with a ruling-line node (labels(:,2)) take that node's folded
% position directly; all others are barycentrically interpolated within
% their enclosing ruling triangle, using the planar (pre-projection) coordinates.
vert_f = zeros(length(vert_p), 3);

for i = 1:length(vert_p)
    if ~isnan(labels(i, 2))
        vert_f(i, 1:3) = ruling_nodes_f(:, labels(i, 2));
    else
        tri_idx = labels(i, 1);
        p = vert_p(i, 1:2);
        p1 = [ruling_nodes_p(1:2, ruling_tris(tri_idx, 1)); 0];
        p2 = [ruling_nodes_p(1:2, ruling_tris(tri_idx, 2)); 0];
        p3 = [ruling_nodes_p(1:2, ruling_tris(tri_idx, 3)); 0];
        p1_f = ruling_nodes_f(1:3, ruling_tris(tri_idx, 1));
        p2_f = ruling_nodes_f(1:3, ruling_tris(tri_idx, 2));
        p3_f = ruling_nodes_f(1:3, ruling_tris(tri_idx, 3));

        dem = (p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*p3(1) + p2(2)*p3(1) + p1(1)*p3(2) - p2(1)*p3(2));
        beta_tri = -(p1(2)*p3(1) - p1(1)*p3(2) - p1(2)*p(1) + p3(2)*p(1) + p1(1)*p(2) - p3(1)*p(2))/dem;
        delta_tri = (p1(2)*p2(1) - p1(1)*p2(2) - p1(2)*p(1) + p2(2)*p(1) + p1(1)*p(2) - p2(1)*p(2))/dem;

        vert_f(i, 1:3) = p1_f + beta_tri*(p2_f-p1_f)+delta_tri*(p3_f-p1_f);
    end
end
end

function alphas = getAlphas(nodes, faces)
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

function tris = makeLadderTris(nCols)
% Two triangles per rung connecting column k of one ruling line (node index
% k) to column k of a second ruling line (node index nCols+k) - a ladder
% strip with no shared apex, unlike the gore ruling fans.
tris = zeros(2*(nCols-1), 3);
for k = 1:(nCols-1)
    tris(2*k-1, :) = [k, k+1, nCols+k];
    tris(2*k, :)   = [k+1, nCols+k+1, nCols+k];
end
end

function vert_p = projectBrimNodes(vert, R, surf_func, brim_func)
% Projects (x,y) onto z for an Nx3 array of brim nodes. Anything exactly at
% radius R is the seam shared with the base mesh (the whole reused
% end_edge/end_edge_gore1 arcs are at r=R, not just the two corners) and
% must use surf_func to match the base mesh bit-for-bit - brim_func(R) is
% not guaranteed to equal surf_func(R) exactly, which the coordinate-based
% exact-match merge in mergeMeshPiece requires. Anything strictly beyond R
% uses brim_func.
r = sqrt(vert(:, 1).^2 + vert(:, 2).^2);
on_seam = abs(r - R) < 1e-9;
vert_p = vert;
vert_p(~on_seam, 3) = brim_func(r(~on_seam));
vert_p(on_seam, 3) = surf_func(r(on_seam));
end

function [vert_u, vert_f, vert_ref, vert_map] = mergeMeshPiece(vert_u, vert_f, vert_ref, vert_p_piece, vert_f_piece, vert_ref_piece)
% Appends a triangulated mesh piece's vertices into the combined mesh,
% reusing any vertex that exactly coincides (by projected x,y,z) with an
% existing one instead of duplicating it. vert_map(i) is the combined-mesh
% index of the piece's local vertex i, for remapping its own face list.
nPiece = size(vert_p_piece, 1);
vert_map = zeros(nPiece, 1);
for i = 1:nPiece
    found_idx = find(ismember(vert_u, vert_p_piece(i, :), 'rows'));
    if isempty(found_idx)
        vert_u = [vert_u; vert_p_piece(i, :)];
        vert_f = [vert_f; vert_f_piece(i, :)];
        vert_ref = [vert_ref; vert_ref_piece(i, :)];
        vert_map(i) = length(vert_u);
    else
        vert_map(i) = found_idx(1);
    end
end
end
