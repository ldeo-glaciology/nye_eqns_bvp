% Run Solver
[s,t,h,Q,S,N,u] = nf_solver();

%% Plotting
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Making plots
figure('Name','Evolution of Drainage');
subplot(3,1,1);
plot(t,h);
xlabel('time, \emph{t} (years)');
ylabel('lake height, \emph{h} (m)');
title('Lake height over time');

subplot(3,1,2);
plot(t,Q(1,:),'DisplayName','Lake Exit');
hold on;
plot(t,Q(end,:),'--','DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)');
ylabel('channel flow rate, \emph{Q} $(m^{3} s^{-1})$');
title('Flow over time');
legend('Location','northwest');

subplot(3,1,3);
plot(t,S(1,:),'DisplayName','Lake Exit'); % now plotting area at lake
hold on;
plot(t,S(end,:),'DisplayName','Channel Exit');
xlabel('time, \emph{t} (years)');
ylabel('channel cross-sectional area, \emph{S} $(m^2)$');
title('Channel area over time');
legend('Location','northwest');

figure('Name', 'Colormaps');
subplot(4,1,1);
imagesc(t,s,N);
xlabel('time, \emph{t} (years)');
ylabel('position along channel, \emph{x} (km)');
title('Effective Pressure Colormap');
cbarN = colorbar;
cbarN.Label.Interpreter = 'latex';
cbarN.TickLabelInterpreter = 'latex';
ylabel(cbarN, "effective pressure, \emph{N} (Pa)")

subplot(4,1,2);
imagesc(t,s,Q);
xlabel('time, \emph{t} (years)');
ylabel('position along channel, \emph{x} (km)');
title('Channel Flow Colormap');
cbarQ = colorbar;
cbarQ.Label.Interpreter = 'latex';
cbarQ.TickLabelInterpreter = 'latex';
ylabel(cbarQ, "channel flow rate, \emph{Q} $(m^3 s^{-1})$")

subplot(4,1,3);
imagesc(t,s,S);
xlabel('time, \emph{t} (years)');
ylabel('position along channel, \emph{x} (km)');
title('Channel Cross-sectional Area Colormap');
cbarS = colorbar;
cbarS.Label.Interpreter = 'latex';
cbarS.TickLabelInterpreter = 'latex';
ylabel(cbarS, "channel cross-sectional area, \emph{S} $(m^2)$")

subplot(4,1,4);
imagesc(t,s,u,[-100,1000]);
xlabel('time, \emph{t} (years)');
ylabel('position along channel, \emph{x} (km)');
title('Glacier Velocity Colormap');
cbaru = colorbar;
cbaru.Label.Interpreter = 'latex';
cbaru.TickLabelInterpreter = 'latex';
ylabel(cbaru, "glacier basal velocity, \emph{u} $(m y^{-1})$");
%set(gca,'ColorScale','log');
