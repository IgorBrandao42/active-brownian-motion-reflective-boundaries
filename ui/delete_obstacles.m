function delete_obstacles()

global user_obstacle_all h_obstacles
user_obstacle_all = user_obstacle_all(1);

for i=2:length(h_obstacles)
    delete(h_obstacles{i});
end

end