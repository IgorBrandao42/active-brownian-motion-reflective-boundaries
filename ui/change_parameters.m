function change_parameters()
global eta T timestamps 

coords = inputdlg({'\eta [N*s/m^2]','T [K]', 'End time [s]', 'Time resolution'}, 'Parameters' , 1, {num2str(eta),num2str(T), num2str(timestamps(end)), num2str(length(timestamps))});

eta      = str2double(coords(1));
T        = str2double(coords(2));
t_end    = str2double(coords(3));
N_time   = str2double(coords(4));

timestamps = linspace(0, t_end, N_time);

% global user_particle_all
% for i=1:N_ensemble
%   x0 = user_particle_all(i).x(1);
%   y0 = user_particle_all(i).y(1);
%   phi0 = user_particle_all(i).phi(1);
%   R0 = user_particle_all(i).R;
%   v0 = user_particle_all(i).v;
%   omega0 = user_particle_all(i).omega;
% 
%   user_particle_all(i) = particle(x0, y0, phi0, R0, v0, omega0, T, eta);  % Update T and eta, but let every other parameter unaltered
% end

end