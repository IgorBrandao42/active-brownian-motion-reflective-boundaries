function plot_obstacles()

global user_obstacle_all

hold on

for i=1:length(user_obstacle_all)
  user_obstacle_all(i).show();
end

end