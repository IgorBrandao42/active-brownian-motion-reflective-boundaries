function plot_trajectories()

global user_particle_all

hold on

for i=1:length(user_particle_all)
  user_particle_all(i).show();
end

end