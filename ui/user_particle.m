function user_particle()

disp("User selection of obstacle!")              % Tell the user what is happening
disp("Click to anywhere inside plot to place a particle")

xy = [];                                         % Variable to store y coords. of the boundary
xy = ginput(1);                                  % Get first point

%% Tenho que deixar o usuário escolher a angulação inicial, raio, velocidade e freq. angular !
user_particle = particle(xy(1), xy(2), 0, 1e-6, 0, 0);
user_particle.show();                            % Show particle

global user_particle_all
user_particle_all = [user_particle_all, user_particle];

end