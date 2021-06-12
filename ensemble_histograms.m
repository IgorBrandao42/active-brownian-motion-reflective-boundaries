function ensemble_histograms()

global user_particle_all ax_main

limits_x = ax_main.XLim;
limits_y = ax_main.YLim;

figure('Name', "Density distribution")
clf

N_ensemble = length(user_particle_all);
N_time = length(user_particle_all(1).x);

x_total = zeros( [N_ensemble, N_time] );
y_total = zeros( [N_ensemble, N_time] );

for i=1:N_ensemble
  x_total(i,:) = user_particle_all(i).x;
  y_total(i,:) = user_particle_all(i).y;
end

x_total_linear = reshape(x_total, N_ensemble*N_time, 1);
y_total_linear = reshape(y_total, N_ensemble*N_time, 1);


%% Plot density plot
subplot(4,4, [5,6,7,9,10,11,13,14,15]);

hist_2d = histogram2(x_total_linear, y_total_linear, 20, 'DisplayStyle','tile','ShowEmptyBins','on');
shading interp
view(0, 90)
% colorbar('southoutside')
xlim(limits_x)
ylim(limits_y)

x = hist_2d.XBinEdges;
x = ( x(2:end) + x(1:end-1) )/2;

y = hist_2d.YBinEdges;
y = ( y(2:end) + y(1:end-1) )/2;

%% Plot histograms of trajectories
subplot(4,4, [1,2,3]);
hist_x = sum(hist_2d.Values, 1);
plot(x, hist_x);
xlim(limits_x)
title("Mean histogram along x direction")
set(gca,'xtick',[])
set(gca,'xticklabel',[])

subplot(4,4, [8,12,16]);
hist_y = sum(hist_2d.Values, 2);
plot(hist_y, y);
ylim(limits_y)
title("Mean histogram along y direction")
set(gca,'ytick',[])
set(gca,'yticklabel',[])


end