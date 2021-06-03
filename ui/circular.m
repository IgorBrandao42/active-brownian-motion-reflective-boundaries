function circular()

global h_obstacles user_obstacle_all

var_T = 8.3e-7;
interior_is_inside = false;

theta = linspace(0, 2*pi, 100);
x_bound_circle = var_T*cos(theta);
y_bound_circle = var_T*sin(theta);

bounding_box = obstacle(x_bound_circle, y_bound_circle, interior_is_inside);

delete(h_obstacles{1})
h_obstacles{1} = bounding_box.show();

if length(user_obstacle_all) == 1
  user_obstacle_all = bounding_box;
else
  user_obstacle_all(1) = bounding_box;
end
end

