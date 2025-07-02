% returns the angle between two 3d vectors

function angle = angleBetweenVectors3d(v1, v2)
    angle = acos(dot(v1, v2)/norm(v1)/norm(v2));
end