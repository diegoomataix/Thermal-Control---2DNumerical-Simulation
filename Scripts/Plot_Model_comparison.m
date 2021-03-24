% clc; clear all;
%% Data

dx=0.14;
T_b = 298.15;
m = 8e0;                                                    % Spatial Subdivisions
M = 14*m;                                                   % Total n of spatial subdivisions
x = linspace(0, dx, M+1);                                     % Spatial Coordinates [m]
%%%%%


%% Comparar E(perfil central) con D
figure()
hold on
myplot(x, SAVED)
myplot1(X, T_stat_plot_central)
axis([0 dx T_b*0.95 max(T_stat_plot_central)*1.05])
ylabel('{\it T} [K]')
xlabel('{\it x} [m]');
legend('Modelo 1D estacionario', 'Modelo 2D estacionario, perfil central')
hold off

%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting function
function myplot(x, y)
plot(x,y, '-k','LineWidth',1)
box on
grid on
grid minor
axis tight
set(gca,'FontSize',18)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting function
function myplot1(x, y)
plot(x,y, '--k','LineWidth',1)
box on
grid on
grid minor
axis tight
set(gca,'FontSize',18)
end