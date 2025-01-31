function PlotCone(Handle_Axis, x, y, z, point_index, log_z, arrow_length, h_radius, v_radius, Color)
    % This function plots a 3D curve with a directional elliptical cone (arrowhead)
    % at the specified location along the curve. It also supports log scaling on the z-axis.
    %
    % Inputs:
    %   x, y, z       - Vectors representing the 3D curve coordinates
    %   point_index   - Index of the point along the curve to place the cone
    %   log_z         - Boolean, set to true for logarithmic scaling on z-axis
    %   h_radius      - Radius along the horizontal direction (x-axis of cone base)
    %   v_radius      - Radius along the vertical direction (y-axis of cone base)
    
    % Apply log scale to z data if specified
    if log_z
        z = log10(z);  % Transform z data to log scale for consistent cone placement
    end

    % Get coordinates for the selected point and the previous point
    x_int = x(point_index); y_int = y(point_index); z_int = z(point_index);
    x_prev = x(point_index-1); y_prev = y(point_index-1); z_prev = z(point_index-1);

    % Calculate direction vector for cone alignment
    direction = -[x_int - x_prev, y_int - y_prev, z_int - z_prev];
    direction = direction / norm(direction);

    % Cone parameters
    n_points = 100;

    % Generate an ellipse in the XY plane for the base of the cone
    theta = linspace(0, 2*pi, n_points);
    ellipse_x = h_radius * cos(theta);   % Horizontal radius
    ellipse_y = v_radius * sin(theta);   % Vertical radius
    cone_x = [zeros(1, n_points); ellipse_x];
    cone_y = [zeros(1, n_points); ellipse_y];
    cone_z = [zeros(1, n_points); arrow_length * ones(1, n_points)];

    % Rotation matrix to align cone with direction vector
    rotation_axis = cross([0, 0, 1], direction);
    rotation_angle = acos(dot([0, 0, 1], direction));
    if norm(rotation_axis) > 0
        R = axang2rotm([rotation_axis / norm(rotation_axis), rotation_angle]);
        rotated_cone = R * [cone_x(:)'; cone_y(:)'; cone_z(:)'];
        cone_x = reshape(rotated_cone(1, :), size(cone_x)) + x_int;
        cone_y = reshape(rotated_cone(2, :), size(cone_y)) + y_int;
        cone_z = reshape(rotated_cone(3, :), size(cone_z)) + z_int;
    else
        % If already aligned, just translate
        cone_x = cone_x + x_int;
        cone_y = cone_y + y_int;
        cone_z = cone_z + z_int;
    end

    % Draw the 3D elliptical cone (arrowhead) with `surf`
    surf(Handle_Axis, cone_x, cone_y, cone_z, 'FaceColor', Color , 'EdgeColor', 'none');
end
