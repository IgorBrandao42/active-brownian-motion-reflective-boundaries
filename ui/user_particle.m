function user_particle()

global h_particles_initial user_particle_all

ButtonName = questdlg('What kind of particle?', 'Particle selection', 'Brownian motion', 'Microswimmer', 'Ensemble', 'Brownian motion');

switch ButtonName
  case 'Microswimmer'
    coords = inputdlg({'\phi_0','R', 'v', '\Omega'}, 'Microswimmer parameters' , 1, {'0','1e-6', '4e-6', '0'});
    
    phi0       = str2double(coords(1));
    R0         = str2double(coords(2));
    v0         = str2double(coords(3));
    omega0     = str2double(coords(4));
    N_ensemble = 1;
    
  case 'Brownian motion'
    phi0       = 0;
    R0         = 1e-6;
    v0         = 0;
    omega0     = 0;
    N_ensemble = 1;
    
  case 'Ensemble'
    coords = inputdlg({'\phi_0','R', 'v', '\Omega', 'N_ensemble'}, 'Microswimmer parameters' , 1, {'0','1e-6', '0', '0', '200'});
    
    phi0       = str2double(coords(1));
    R0         = str2double(coords(2));
    v0         = str2double(coords(3));
    omega0     = str2double(coords(4));
    N_ensemble = str2double(coords(5));
  otherwise
    disp('Unknown method for particle creation!')
    return
end

disp("User selection of obstacle!")              % Tell the user what is happening
disp("Click to anywhere inside plot to place a particle")

% xy = [];                                         % Variable to store y coords. of the boundary
xy = ginput(1);                                  % Get first point


for i=1:N_ensemble
  user_particle = particle(xy(1), xy(2), phi0, R0, v0, omega0);                 % Create instance of particle
  h_particles_initial{end+1} = user_particle.show();                            % Show particle
  
  user_particle_all = [user_particle_all, user_particle];
end

end


