% I am using intersections: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
% Remember to cite!

fig = figure(2)
clf

x0     = -0.1e-6;                                     % Initial x coordinate [m]
y0     = +0.0e-6;                                     % Initial y coordinate [m]
phi0   = 0;                                      % Initial orientation[rad]
R0     = 1e-6;                                   % Particle radius[m]
v0     = 1e-6;                                   % Self-propulsion velocity [m/s]
omega0 = 0;                                     % [rad/s] (omega>0 -> anti-clockwise; omega<0 -> clockwise)

%color0 = [0,0,1];

subject = particle(x0, y0, phi0, R0, v0, omega0);

t = linspace(0, 6.25/2, 1e4);

var_T = sqrt( subject.D_T*(t(end)-t(1)) );
% Rectangle
x_bound = [-var_T, -var_T, +var_T, +var_T, -var_T];
y_bound = [-var_T, +var_T, +var_T, -var_T, -var_T];

% var_T = 0.1*var_T;
% x_bound_interior = [-var_T, -var_T, +var_T, +var_T, -var_T];
% y_bound_interior = [-var_T, +var_T, +var_T, -var_T, -var_T];

% Triangle
% x_bound = [-var_T, -var_T, +var_T, -var_T];
% y_bound = [-var_T, +var_T, +var_T, -var_T];

interior_is_inside = false;
bound = obstacle(x_bound, y_bound, interior_is_inside);

clf
bound.show()
hold on
plot(subject.x(1), subject.y(1), 'k', 'Marker', '*')

% [x_bound_interior, y_bound_interior] = ginputc('Color', 'k');
% mousePointCoords = ginput();
% x_bound_interior = mousePointCoords(:,1);
% y_bound_interior = mousePointCoords(:,2);
% 
% bound_interior = obstacle(x_bound_interior, y_bound_interior, true);
% bound_interior.show()

% Check if the boundary and unit normal vectors are correct !
% figure(1)
% plot(subject.x(1), subject.y(1), 'k', 'Marker', '*')
% hold on
% bound.show()
% hold off

subject.time_evolution(t, [bound, bound_interior])


plot(subject.x, subject.y)
plot(subject.x(1), subject.y(1), 'k', 'Marker', '*')
plot(subject.x(end), subject.y(end), 'r', 'Marker', '*')


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


% end




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



