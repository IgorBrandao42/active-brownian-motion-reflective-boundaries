% I am using intersections: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections

%% Define plotting figure
fig = figure(1);
fig.Units = 'pixels';
fig.Position = [599.4000  180.2000  559.2000  502.4000];
                 % subtightplot(m, n, p, gap        , marg_h    , marg_w    ,varargin)
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.1 0.01], [0.1 0.01]);
clf
% subplot(4,4, [5,6,7,9,10,11,13,14,15]);

%% Parameters for microswimmer and timestamps for simulation
x0     = -0.0e-6;                                % Initial x coordinate [m]
y0     = +0.0e-6;                                % Initial y coordinate [m]
phi0   = 0;                                      % Initial orientation[rad]
R0     = 5e-6;                                   % Particle radius[m]
v0     = 1e-6;                                   % Self-propulsion velocity [m/s]
omega0 = 0;                                      % [rad/s] (omega>0 -> anti-clockwise; omega<0 -> clockwise)
% color0 = [0,0,1];                              % Color for the particle (if this parameters is not given to the particle contructor, a random colo will be assigned!)

t = linspace(0, 10*6.25, 3e3);                      % Timestampes for simulations


%% Fluid parameters
T   = 300;                                 % Environmental temperature [K]
eta = 0.001;                               % Fluid viscosity           [N*s/m^2]

%% Define the particle
subject = [];
subject = [subject, particle(x0, y0, phi0, R0, v0, omega0, T, eta, [0.8,0.8,0])];  % Create particle
subject(1).show();                                    % Show initial position of the particle
hold on

subject = [subject, particle(x0, y0, phi0, R0, v0, omega0, T, eta, [0,1,0])];  % Create particle
subject(2).show();                                    % Show initial position of the particle

subject = [subject, particle(x0, y0, phi0, R0, v0, omega0, T, eta, [0,0,1])];  % Create particle
subject(3).show();                                    % Show initial position of the particle


%% Rectangular bounding box
var_T = 5e-6; %sqrt( subject.D_T*(t(end)-t(1)) );
x_bound = [-var_T, -var_T, +var_T, +var_T, -var_T];
y_bound = [-var_T, +var_T, +var_T, -var_T, -var_T];
interior_is_inside = false;

bound = obstacle(x_bound, y_bound, interior_is_inside); % Create bounding box
bound.show()                                       % Plot bounding box
hold on


%% Inner rectangular obstacle
var_T_inner = 0.1*var_T;
% x_bound_interior = [-var_T_inner, -var_T_inner, +var_T_inner, +var_T_inner, -var_T_inner];
% y_bound_interior = [-var_T_inner, +var_T_inner, +var_T_inner, -var_T_inner, -var_T_inner];

% x_bound_interior = -var_T_inner; % , -var_T_inner, +var_T_inner, +var_T_inner, -var_T_inner];
% y_bound_interior = -var_T_inner; % , +var_T_inner, +var_T_inner, -var_T_inner, -var_T_inner]
% 
theta = linspace(0, pi, 50);
x_bound_interior = var_T_inner + 6*var_T_inner*sin(theta);
y_bound_interior = 3*var_T_inner*cos(theta);
% 
% x_bound_interior(end+1) = x_bound_interior(1);
% y_bound_interior(end+1) = y_bound_interior(1);
% 
x_bound_interior = [x_bound_interior, 4*var_T_inner*sin(pi-theta)];
y_bound_interior = [y_bound_interior, 3*var_T_inner*cos(pi-theta)];

bound_interior = obstacle(x_bound_interior, y_bound_interior, true); % Create inner obstacle
bound_interior.show()                              % Plot obstacle


for i=1:length(subject)
  subject(i).time_evolution(t, [bound, bound_interior])
  subject(i).show();
end

% plot(subject.x, subject.y)                         % Plot trajectories
xlim([-var_T, +var_T])
ylim([-var_T, +var_T])

xlabel('x [\mum]')
ylabel('y [\mum]')


%% Adjust labels
ax = gca;
ax.FontSize = 25;

xticklabels = round(1e6*ax.XTick, 1);
ax.XTickLabel = xticklabels;

yticklabels = round(1e6*ax.YTick, 1);
ax.YTickLabel = yticklabels;

xlabel('x [\mum]');
ylabel('y [\mum]');


%% Plot histograms of trajectories
% subplot(4,4, [1,2,3]);
% hist_x = histogram(subject.x);
% xlim([-var_T, +var_T])
% title("Histogram along x direction")
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% 
% subplot(4,4, [8,12,16]);
% hist_y = histogram(subject.y, 'Orientation', 'horizontal');
% ylim([-var_T, +var_T])
% title("Histogram along y direction")
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])

% 
% A = [0 0 1 1 1 0 0 0 0 NaN NaN 1 0 0 0 1 0 1 0 1 0 0 0 1 1 1 1];
% C = categorical(A,[1 0 NaN],{'yes','no','undecided'})
