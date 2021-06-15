function load_canvas()

global user_obstacle_all user_particle_all h_obstacles h_particles_initial h_trajectories

[baseFileName, folderName, ~] = uigetfile('*.mat');
if ~isequal(baseFileName, 0)
  inputFullFileName = fullfile(folderName, baseFileName);
else
  return;
end

delete(h_obstacles{1});

load(inputFullFileName)

h_obstacles{1} = user_obstacle_all(1).show();


plot_obstacles()
plot_trajectories()

end




% Get file_name MATLAB UI, see RPN Philio
% output_directory = 'canvas_configuration';
% if ~exist(output_directory, 'dir')
%   mkdir(output_directory);
% end

% input_directory = 'Video_Input';
% addpath(input_directory);
% addpath(output_directory);                       % Add path to output's folder



