function [ellipsoidalCoords] = CartesianToEllipsoid(cartesianCoords, center, algebraic)
    x = cartesianCoords(:, 1);
    y = cartesianCoords(:, 2);
    z = cartesianCoords(:, 3);

    x_shifted = x - center(1);
    y_shifted = y - center(2);
    z_shifted = z - center(3);

    
end

