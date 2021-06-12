function simulate()

global user_obstacle_all user_particle_all

t = linspace(0, 6.25/2, 1e3);

for i=1:length(user_particle_all)
  
  user_particle_all(i).time_evolution(t, user_obstacle_all)
  
end


end