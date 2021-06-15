function simulate()

global user_obstacle_all user_particle_all timestamps

for i=1:length(user_particle_all)
  user_particle_all(i).time_evolution(timestamps, user_obstacle_all)
end


end