function ensemble_histograms()

global user_particle_all ax_main

ax = ax_main;
h_main_plot = subplot(4,4, [5,6,7,9,10,11,13,14,15], ax); % Transforms main axis into subplot !
limits_x = h_main_plot.XLim;
limits_y = h_main_plot.YLim;

x_mean = zeros(size(user_particle_all(1).x));
y_mean = zeros(size(user_particle_all(1).y));

N_ensemble = length(user_particle_all);
for i=1:N_ensemble
  x_mean = x_mean + user_particle_all(i).x;
  y_mean = y_mean + user_particle_all(i).y;
end

x_mean = x_mean/N_ensemble;
y_mean = y_mean/N_ensemble;

%% Plot histograms of trajectories
subplot(4,4, [1,2,3]);
histogram(x_mean);
xlim(limits_x)
title("Histogram along x direction")
set(gca,'xtick',[])
set(gca,'xticklabel',[])

subplot(4,4, [8,12,16]);
histogram(y_mean, 'Orientation', 'horizontal');
ylim(limits_y)
title("Histogram along y direction")
set(gca,'ytick',[])
set(gca,'yticklabel',[])

ax_main = h_main_plot;


end