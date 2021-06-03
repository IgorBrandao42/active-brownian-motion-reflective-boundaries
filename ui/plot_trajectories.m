function plot_trajectories()

global user_particle_all h_trajectories

hold on

for i=1:length(user_particle_all)
  h_trajectories{i} = user_particle_all(i).show();
end

end