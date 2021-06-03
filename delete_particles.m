function delete_particles()

global user_particle_all h_trajectories
user_particle_all = [];

for i=1:length(h_trajectories)
    delete(h_trajectories{i});
end

end