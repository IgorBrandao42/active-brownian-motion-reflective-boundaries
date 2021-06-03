function clear_bounding_box

global h_trajectories
for i=1:length(h_trajectories)
    delete(h_trajectories{i});
end

global h_obstacles

for i=2:length(h_obstacles)
    delete(h_obstacles{i});
end

end