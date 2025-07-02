% returns the 2D counterclockwise angle in radians between vector1 and vector2
% (rotating vector1 counterclockwise by angle will give vector 2)
% the output is between -pi and pi

function angle = angleBetweenVectors(vector1, vector2)
    angle = atan2(vector2(2), vector2(1)) - atan2(vector1(2), vector1(1));
    
    if angle > pi()
        angle = angle - 2*pi();
    end
    
    if angle < -pi()
        angle = angle + 2*pi();
    end
end