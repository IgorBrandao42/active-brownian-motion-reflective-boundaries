function rectangular()

global h_obstacles user_obstacle_all

var_T = 8.3e-7;
interior_is_inside = false;

x_bound = [-var_T, -var_T, +var_T, +var_T, -var_T];
y_bound = [-var_T, +var_T, +var_T, -var_T, -var_T];
bounding_box = obstacle(x_bound, y_bound, interior_is_inside);

delete(h_obstacles{1});
h_obstacles{1} = bounding_box.show();

if length(user_obstacle_all) == 1
  user_obstacle_all = bounding_box;
else
  user_obstacle_all(1) = bounding_box;
end

end

