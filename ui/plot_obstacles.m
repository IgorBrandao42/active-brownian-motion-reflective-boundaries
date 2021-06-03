function plot_obstacles()

global user_obstacle_all h_obstacles

hold on
h_obstacles = cell(length(user_obstacle_all), 1);
for i=2:length(user_obstacle_all)
  h_obstacles{i} = user_obstacle_all(i).show();
end

end