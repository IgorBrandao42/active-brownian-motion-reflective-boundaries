function user_obstacle()

global fig

disp("User selection of obstacle!")              % Tell the user what is happening
disp("Click to draw obstacle boundary and press key 'o' to stop")

max_points = 1e2;                                % Maximum number of points the user can choose (I don't want an infinite loop!) 

xy = zeros(2, max_points);                       % Variable to store y coords. of the boundary
xy(:, 1) = ginput(1);                            % Get first point
h0 = plot(xy(1, 1), xy(2, 1), 'r', 'Marker', '*'); % Plot first point

h1 = [];
for i=2:max_points                               % Loop until user presses key and stops drawing
  
  temp = ginput(1);                              % Wait user to click next point and get its coords
  
  if get(fig,'CurrentCharacter') == 'o'          % If a key was pressed ( https://www.mathworks.com/matlabcentral/answers/57157-how-can-i-get-matlab-to-leave-a-loop-when-any-key-is-pressed )
    break                                        % End loop
  end
  
  xy(:, i) = temp;                               % Workaround ginput not returning a 2-vector because user stopped the loop !
  
  delete(h1)
  h1 = plot(xy(1, 1:i), xy(2, 1:i), 'r');        % Plot every point choosen (the boundary so far)
  drawnow
  
end
delete(h0)                                       % Delete first point selected
delete(h1)                                       % Delete boundary so far

fig.CurrentCharacter = 'q';

user_obstacle = obstacle(xy(1, 1:i), xy(2, 1:i), true); % Create obstacle
user_obstacle.show();                                % Show obstacle

global user_obstacle_all
user_obstacle_all = [user_obstacle_all, user_obstacle];

end