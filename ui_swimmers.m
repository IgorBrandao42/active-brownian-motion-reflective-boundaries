%% TO DO LIST
% Remember to cite intersections: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections

% Create custom particles (different radii, self propulsion velocities, ...)
% Custom plotting region
% Custom timestamps

%% Set up environment
addpath('ui')

global user_obstacle_all user_particle_all fig h_trajectories h_obstacles h_particles_initial ax_main timestamps eta T
user_obstacle_all = [];
user_particle_all = [];
h_trajectories = {};
h_obstacles = {};
h_particles_initial = {};

timestamps = linspace(0, 6.25, 1e3);       % Timestamps for the simulation [s]
T   = 300;                                 % Environmental temperature     [K]
eta = 0.001;                               % Fluid viscosity               [N*s/m^2]

%% Create figure
fig = figure('Name', "Microswimmer simulation");
clf
% fig.Color = [0 0.5 0.5];

fig.MenuBar = 'none';                           % Hide default menubar
fig.ToolBar = 'figure';                         % Show toolbar


%% Create menus
menu_add = uimenu('Label', 'Add to canvas');
menu_add_obstacle = uimenu(menu_add, 'Label', 'Add obstacle', 'Callback', 'user_obstacle');
menu_add_particle = uimenu(menu_add, 'Label', 'Add particle', 'Callback', 'user_particle');

menu_simulation = uimenu('Label', 'Simulate');
menu_simulate   = uimenu(menu_simulation, 'Label', 'Simulate', 'Callback', 'simulate');

menu_show = uimenu('Label', 'Show');
menu_show_obstacles = uimenu(menu_show, 'Label', 'Show obstacles'   , 'Callback', 'plot_obstacles');
menu_show_traject   = uimenu(menu_show, 'Label', 'Show trajectories', 'Callback', 'plot_trajectories');
menu_animate        = uimenu(menu_show, 'Label', 'Show animation'   , 'Callback', 'animate_particles');
menu_histograms     = uimenu(menu_show, 'Label', 'Show ensemble histograms', 'Callback', 'ensemble_histograms');

menu_clear = uimenu('Label','Clear');
menu_clear_canvas     = uimenu(menu_clear, 'Label', 'Clear canvas'    , 'Callback', 'clear_bounding_box');
menu_delete_particles = uimenu(menu_clear, 'Label', 'Delete particles', 'Callback', 'delete_particles');
menu_delete_obstacles = uimenu(menu_clear, 'Label', 'Delete obstacles', 'Callback', 'delete_obstacles');

menu_save_load = uimenu('Label','Save/load');
menu_load_canvas      = uimenu(menu_save_load, 'Label', 'Load obstacles and particles', 'Callback', 'load_canvas');
menu_save_canvas      = uimenu(menu_save_load, 'Label', 'Save obstacles and particles', 'Callback', 'save_canvas');

menu_boundary  = uimenu('Label', 'Change bounding box');
menu_rectangular_boundary = uimenu(menu_boundary, 'Label', 'Retangular', 'Callback', 'rectangular');
menu_circular_boundary    = uimenu(menu_boundary, 'Label', 'Circular'  , 'Callback', 'circular');

menu_environment  = uimenu('Label', 'Change parameters');
menu_env          = uimenu(menu_environment, 'Label', 'Change parameters', 'Callback', 'change_parameters');

%% Plotting region (bounding box)
var_T = 8.3e-7;
interior_is_inside = false;

x_bound = [-var_T, -var_T, +var_T, +var_T, -var_T];
y_bound = [-var_T, +var_T, +var_T, -var_T, -var_T];
bounding_box = obstacle(x_bound, y_bound, interior_is_inside);

h_obstacles{1} = bounding_box.show();
hold on

xlim(1.01*[-var_T, +var_T]);
ylim(1.01*[-var_T, +var_T]);
user_obstacle_all = bounding_box;

ax_main = gca;

axis square






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


% user_obst = user_obstacle(fig);

% plot(subject.x, subject.y)
% plot(subject.x(end), subject.y(end), 'r', 'Marker', '*')

% x0     = -0.1e-6;                                     % Initial x coordinate [m]
% y0     = +0.0e-6;                                     % Initial y coordinate [m]
% phi0   = 0;                                      % Initial orientation[rad]
% R0     = 1e-6;                                   % Particle radius[m]
% v0     = 1e-6;                                   % Self-propulsion velocity [m/s]
% omega0 = 0;                                     % [rad/s] (omega>0 -> anti-clockwise; omega<0 -> clockwise)
%color0 = [0,0,1];
% subject = particle(x0, y0, phi0, R0, v0, omega0);
% subject.show();

% Small square
% var_T = 0.1*var_T;
% x_bound_interior = [-var_T, -var_T, +var_T, +var_T, -var_T];
% y_bound_interior = [-var_T, +var_T, +var_T, -var_T, -var_T];

% Triangle
% x_bound = [-var_T, -var_T, +var_T, -var_T];
% y_bound = [-var_T, +var_T, +var_T, -var_T];

% N_time = 1e3;
% t = linspace(0, 6.25, N_time);
% 
% dt = t(2) - t(1);
% 
% 
% 
% 
% w_x   = randn([1, N_time]);
% w_y   = randn([1, N_time]);
% w_phi  
% A_R = sqrt(2*subject.D_R*dt);
% A_T = sqrt(2*subject.D_T*dt);
% 
% subject.x    = [subject.x(1)  , zeros([1, N_time-1])];         % Variable to store x position  at each time
% subject.y    = [subject.y(1)  , zeros([1, N_time-1])];         % Variable to store x position  at each time
% subject.phi  = [subject.phi(1), zeros([1, N_time-1])];         % Variable to store orientation at each time
% 
% for k = 1:N_time-1                             % Euler-Maruyama method for stochastic integration
%   subject.x(k+1)   = subject.x(k) + subject.v*cos(subject.phi(k))*dt + A_T*w_x(k);
%   subject.y(k+1)   = subject.y(k) + subject.v*sin(subject.phi(k))*dt + A_T*w_y(k);
%   
%   subject.phi(k+1) = subject.x(k) + subject.omega*dt + A_R*w_phi(k);
% end  % loop on time




%%%%%%%%%%% Plotting function %%%%%%%%%

% function plot_semi_classical(obj)
% mode_colors = [[0, 0, 0] ; [0, 0, 1]; [1, 0, 0]; [0, 0.5, 0]; [0.75, 0, 0.75]; [0.75, 0.75, 0]; [0.6350, 0.0780, 0.1840]]; % Array with colors to be consistently used throughout the plots
% 
% if isempty(obj.R_semi_classical)
%   disp("No semi-classical quadratures calulated, please calculate it before plotting it!")
%   return
% end
% 
% N_particles = obj.Size_matrices/2.0 - 1.0;
% 
% figure('Name', "Semi-classical Quadratures Phase Space");           % Open a figure to contain the subplots
% 
% x_label = "x_" + string(1:N_particles);          % Create the labels for the x axis
% x_label = ["q", x_label];
% 
% y_label = "p_" + string(1:N_particles);          % Create the labels for the y axis
% y_label = ["p", y_label];
% 
% title_name = "Particle " + string(1:N_particles);% Create the title (each mode name)
% title_name = ["Cavity", title_name];
% 
% for i=0:N_particles                              % For each mode
%   x = obj.R_semi_classical(2*i+1, :);            % Position quadrature for the i-th mode (semi-classical)
%   p = obj.R_semi_classical(2*i+2, :);            % Momentum quadrature for the i-th mode (semi-classical)
%   
%   dq = x(2:end)-x(1:end-1);                          % Increments in the position for each time step to be used in the arrow plot
%   dp = p(2:end)-p(1:end-1);                          % Increments in the momentum for each time step to be used in the arrow plot
%   
%   quiver_scale = 0;                                  % Tell MATLAB to not scale the arrows
%   
%   subplot(1, N_particles+1, i+1);                % Specify why subplot to use
%   quiver(x(1:end-1), p(1:end-1), dq, dp, quiver_scale, 'Color', mode_colors(i+1,:)) % Arrow plot the phase space trajectory
%   hold on
%   plot(x(1), p(1), 'x', 'Color', mode_colors(i+1,:), 'MarkerSize', 12)
%   
%   title(title_name(i+1))                             % Add a title
%   xlabel(x_label(i+1));                              % Label the x axis
%   ylabel(y_label(i+1),'rotation',0,'VerticalAlignment','middle'); % Label the p axis
%   
%   set(gca, 'Fontsize', 18)                           % Set text to be on proper font size
%   
% end
% 
% end



