function [angle, n, m] = dihedralAngle(p_j,p_k,p_i,p_l)
% Vectorized signed dihedral, robust to inconsistent face ordering.
% Inputs are n×3, output theta is n×1.

    % Edge direction (normalize row-wise)
    rij = p_i-p_j;
    rkj = p_k-p_j;
    rkl = p_k-p_l;

    % Raw normals
    m = cross(rij, rkj, 2);
    n = cross(rkj, rkl, 2);

    angle = angleBetweenVectors3d(m, n);
    angle = real(angle);

    eta = ones(size(angle)).*sign(sum(m.*rkl,2));
    eta(sum(m.*rkl,2)==0) = 1;
    angle = eta.*angle;

    angle = mod(angle, 2*pi);

    m = m./vecnorm(m,2,2);
    n = n./vecnorm(n,2,2);

    % % check direction
    % ref = cross(n1, n2, 2); ref = ref ./ vecnorm(ref, 2, 2);
    % angle(dot(e, ref, 2) > 0) = -angle(dot(e, ref, 2) > 0);

    % Compute signed dihedral
    % x = dot(n1, n2, 2);
    % y = dot(e, cross(n1, n2, 2), 2);
    % angle = atan2(y, x);   % n×1
end