function animate_particles()

global user_particle_all fig h_trajectories

disp("Starting animation!")
disp("Click to draw obstacle boundary and press key 'a' to stop")

t = linspace(0, 6.25/2, 1e4);

hold on
h_trajectories = cell(length(user_particle_all),1);

for k=1:length(t)
  for i=1:length(user_particle_all)
    delete(h_trajectories{i});
    h_trajectories{i} = user_particle_all(i).show(k);
    drawnow
  end
  
  if get(fig,'CurrentCharacter') == 'a'       % If a key was pressed ( https://www.mathworks.com/matlabcentral/answers/57157-how-can-i-get-matlab-to-leave-a-loop-when-any-key-is-pressed )
    break                                        % End loop
  end
end


end