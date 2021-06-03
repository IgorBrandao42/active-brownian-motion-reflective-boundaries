% I am using intersections: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections

fig = figure(2)
clf

x0     = -0.1e-6;                                % Initial x coordinate [m]
y0     = +0.0e-6;                                % Initial y coordinate [m]
phi0   = 0;                                      % Initial orientation[rad]
R0     = 1e-6;                                   % Particle radius[m]
v0     = 1e-6;                                   % Self-propulsion velocity [m/s]
omega0 = 0;                                      % [rad/s] (omega>0 -> anti-clockwise; omega<0 -> clockwise)
% color0 = [0,0,1];


% Rectangular bounding box
var_T = sqrt( subject.D_T*(t(end)-t(1)) );
x_bound = [-var_T, -var_T, +var_T, +var_T, -var_T];
y_bound = [-var_T, +var_T, +var_T, -var_T, -var_T];
interior_is_inside = false;

bound = obstacle(x_bound, y_bound, interior_is_inside); % Create bounding box
bound.show()                                       % Plot bounding box

hold on

% Inner rectangular obstacle
var_T = 0.1*var_T;
x_bound_interior = [-var_T, -var_T, +var_T, +var_T, -var_T];
y_bound_interior = [-var_T, +var_T, +var_T, -var_T, -var_T];

bound_interior = obstacle(x_bound_interior, y_bound_interior, true); % Create inner obstacle
bound_interior.show()                              % Plot obstacle

subject = particle(x0, y0, phi0, R0, v0, omega0);  % Create particle
subject.show();                                    % Show initial position of the particle

t = linspace(0, 6.25/2, 1e4);                      % Timestampes for simulations
subject.time_evolution(t, [bound, bound_interior]) % Simulation

plot(subject.x, subject.y)                         % Plot trajectories

