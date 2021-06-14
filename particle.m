classdef particle < handle                 % Class simulating an experiment with N nanoparticles linearly interacting with a single optical cavity field
  properties
    x                                            % x position  at each time
    y                                            % y position  at each time
    phi                                          % orientation at each time
    
    R                                            % Particle radius
    v                                            % Self-propulsion velocity
    omega                                        % Orientation angular velocity
    
    D_T                                          % Translational diffusion coefficient
    D_R                                          % Rotational    diffusion coefficient
    
    color                                        % Color of the particle                (just for the plotting, does not influence anything else!)
  end
  
  methods
    function obj = particle(x0, y0, phi0, R0, v0, omega0, T, eta, color0)
      % Class constructor for a spherical particle with radius R and color color0
      % initially at position (x0,y0) and orientation phi0
      
      obj.x = x0;
      obj.y = y0;
      obj.phi = phi0;
      
      obj.R = R0;
      obj.v = v0;
      obj.omega = omega0;
      
      if nargin < 9
        obj.color = rand(3,1);
      else
        obj.color = color0;
      end
      
      if nargin < 7
        T   = 300;                                 % Environmental temperature [K]
        eta = 0.001;                               % Fluid viscosity           [N*s/m^2]
      end
      
      % Contants for the calculation of the diffusion coefficients
      k_B = 1.380649e-23;                        % Boltzmann's constant      [J/K]
      
      obj.D_T = k_B*T/(6*pi*eta*R0);             % [metre^2/s] = 1e-12*[micrometer^2/s]
      obj.D_R = k_B*T/(8*pi*eta*R0^3);           % [rad^2/s]
    end
    
    function time_evolution(obj, t, obstacles)
      % Solve the overdamped Langevin equations for translational and rotational dynamics
      % of a circular particle in 2D given some obstacles, which the particle can bounce from
      %
      % obj       - instance of particle class
      % t         - timestamps for the simulation
      % osbtacles - array of function handles. GIven some position, returns a boolean telling if the particle is inside the boundary
      %
      %
      N_time = length(t);
      
      dt = t(2) - t(1);
      
      w_x   = randn([1, N_time]);
      w_y   = randn([1, N_time]);
      w_phi = randn([1, N_time]);
      
      A_R = sqrt(2*obj.D_R*dt);
      A_T = sqrt(2*obj.D_T*dt);
      
      obj.x    = [obj.x(1)  , zeros([1, N_time-1])];  % Variable to store x position  at each time
      obj.y    = [obj.y(1)  , zeros([1, N_time-1])];  % Variable to store x position  at each time
      obj.phi  = [obj.phi(1), zeros([1, N_time-1])];  % Variable to store orientation at each time
      
      for k = 1:N_time-1
        obj.x(k+1)   = obj.x(k)   + obj.v*cos(obj.phi(k))*dt + A_T*w_x(k);  % Euler-Maruyama method for stochastic integration
        obj.y(k+1)   = obj.y(k)   + obj.v*sin(obj.phi(k))*dt + A_T*w_y(k);
        obj.phi(k+1) = obj.phi(k) + obj.omega*dt             + A_R*w_phi(k);
        
        obj.handle_borders(obstacles, k);             % Test if the new position is within an obstacle and correct the last position (reflection)
      end
    end
    
    function handle_borders(obj, obstacles, k)
      reset_loop = true;                              % Make the while loop run at least one time !
        while reset_loop                              % Make sure the particle was not previously reflected into another obstacle !
          reset_loop = false;                         % Only redo the while loop if the particle was reflected
          
          for j=1:length(obstacles)                   % For each boundary
            if obstacles(j).was_penetrated(obj.x(k+1), obj.y(k+1))  % If the particle got inside the boundary
              [obj.x(k+1), obj.y(k+1)] = obstacles(j).reflect(obj, k);           % Reflect the particle on the boundary
              reset_loop = true;                      % Redo everything to test if there it was reflected into another obstacle or in another region of the same obstacle
              break                                   % Break out of inner loop
            end  % if statement
          end  % loop on boundaries
          
          % Tem um problema aqui, numa eventual segunda iteração do while loop, deveria considerar reflexão em relação a interseção do for-loop anterior com a fronteira
          % Mas fazer isso vai ferrar a notação e é um caso muito específico pra um projeto desses !
          
        end  % loop on redoing boundaries
    end
    
    function h = show(obj, end_time)
      
      if nargin == 1
        end_time = length(obj.x);
      end
      
      if length(obj.x) == 1
        h = plot(obj.x(1), obj.y(1), 'Color', obj.color, 'Marker', '*', 'MarkerSize', 10);
      else
        h = plot(obj.x(1:end_time), obj.y(1:end_time), 'Color', obj.color, 'Linewidth', 1);
        
      end
    end
    
  end % End methods
end % End classdef