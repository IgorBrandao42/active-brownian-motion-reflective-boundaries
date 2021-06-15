function save_canvas()

global user_obstacle_all user_particle_all h_obstacles h_particles_initial h_trajectories eta T timestamps

save("canvas_" + datestr(now,'dd-mmmm-yyyy_HH-MM-SS'), "user_obstacle_all", "user_particle_all", "h_obstacles", "h_trajectories", "h_particles_initial", "eta", "T", "timestamps");

end