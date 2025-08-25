% returns the angle between two 3d vectors

function angle = angleBetweenVectors3d(v1, v2)
    % v1 and v2 are n x 3 matrices
    angle = acos(sum(v1.*v2,2)./(vecnorm(v1,2,2) .* vecnorm(v2,2,2)));
end