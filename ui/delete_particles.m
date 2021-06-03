function delete_particles()

global user_particle_all h_trajectories h_particles_initial
user_particle_all = [];

for i=1:length(h_trajectories)
    delete(h_trajectories{i});
end

for i=1:length(h_particles_initial)
    delete(h_particles_initial{i});
end

end