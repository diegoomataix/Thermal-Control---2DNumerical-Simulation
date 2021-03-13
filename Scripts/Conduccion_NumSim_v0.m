%% **********************************************************************************
%                      CONDUCCIÓN DE CALOR, SIMULACIÓN NUMÉRICA
% ------------------------------------------------------------------------------------
% Realizado por Diego Mataix Caballero.
%
%  ADDITIONAL NOTES:
% PCB de FR-4 =: 140 x 100 x 1.5 (dx * dy * dz)
% Recubrimiento de Cu de 50e-6 m
%       - en cara 1 : continuo
%       - en cara 2 : 90% FR-4, 10% Cu
% 3 IC, cada uno disipa 5W, con k_ic = 5 [W/(mK)], con c_ic = 20 [J/K]
%       - distribuidos uniformemente en la PCB, 20 mm de separacion
% PCB tiene contact termico perfecto con paredes permanentemente a 25C, los
% otros dos bordes están térmicamente aislados.
%___________________________________________________________________________
close all; clear all; clc;
%___________________________________________________________________________
%% Datos
Conduccion_NumSim_DATOS
%___________________________________________________________________________
%% Choose exercise to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choose = 'b';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%___________________________________________________________________________
%% Define global parameters
%%% Define coefficients and some parameters
phi = (3 * Q_ic) / Vol;                             % Volumetric dissipation [W/m^3]
l =      [t_rec dz_pcb t_rec];                      % Dimension Vector [m]
k_vect = [k_Cu k_plano (0.1*k_Cu+0.9*k_plano)];     % Conductivity Vector [W/(m·K)] tercera capa es donde van los IC, cubierta solo al 10% de cobre
k_eff = effective(k_vect, l);                       % Effective Conductivity [W/(m·K)]
L = dx/2;                                           % [m]
A = dy * dz;                                        % [m^2]
%%% Mesh
m = 1e2;                                            % Spatial Subdivisions
M = 14*m;                                           % Total n of spatial subdivisions
x = linspace(0, dx, M);                             % Spatial Coordinates [m]
%___________________________________________________________________________
switch(choose)
%___________________________________________________________________________
    case 'a'
        %% Apartado A
        % Considerando que la tarjeta sólo evacua calor por los bordes,
        % determinar la temperatura máxima que se alcanzaría si toda la disipación
        % estuviese uniformemente repartida en la PCB y los IC no influyeran.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        method = 2;         % 1: Only max Temp         2: Show all Temp profile
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch(method)
            case 1
                DT = 1/8 * ( phi * dx^2 / k_eff );                  % Delta T [K]
                T_0 = T_b + DT                                      % Max T [K]
                T_0_C = convtemp(T_0, 'K', 'C')                     % Max T [C]
            case 2
                L = dx/2;                                           % [m]
                b = b_coef(Q_ic_tot, k_eff, dz*dy);                 % [K*m]
                for i = 1:M
                     T(i) = temp_parb(T_b, b, (x(i)), phi, k_eff);  % [K]
                end
                T_0 = max(T)                                        % Max T [K]
                T_0_C = convtemp(T_0, 'K', 'C')                     % Max T [Celsius]
                
                figure()
                myplot(x,T)
                hold on
                axis([0 dx T_b*0.9 max(T)*1.1])
                ylabel('{\it T} [K]')
                xlabel('{\it x} [m]');
                xline(L, '-.')
                yline(T_b, '--')
                yline(T_0, '--')
                hold off
        end
%___________________________________________________________________________
        %% Apartado B
    case 'b'
        % Considerando que la tarjeta sólo evacua calor por los bordes, determinar
        % la temperatura máxima que se alcanzaría con un modelo unidimensional en el
        % que los IC llegaran hasta los bordes aislados, en el límite kIC→∞, y con la kIC dada.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        method = 1;         % 1: Method @ k -->  k_ic         2: Method @ k --> inf
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Límites de cada tramo
        for i = 1:6
            lim(i) = 2*i*m;                                 % [m]
        end
        
        % Define the heat for each type of section: 1: w/o IC, 2 w/ IC
        Q1 = Q_ic_tot;                                      % [W]
        Q2 = Q_ic_tot*(1/3);                                % [W]
        
        % Coefficients for the first section
        a = T_b;                                            % [K]
        b = b_coef(Q_ic_tot, k_eff, dz*dy);                 % [K*m]
        
        % Define the parameters for the sections containing the IC
        l_ic =      [t_rec dz_pcb t_rec dz_ic];                     % Dimension Vector [m]
        k_vect_ic = [k_Cu k_plano (0.1*k_Cu+0.9*k_plano) k_ic];     % Conductivity Vector [W/(m·K)] tercera capa es donde van los IC, cubierta solo al 10% de cobre
        k_eff_ic = effective(k_vect_ic, l_ic);                      % Effective Conductivity [W/(m·K)]
        A_ic = (dz+dz_ic)*dy;                                       % [m^2]                                                        
        phi_ic = Q_ic / (0.02 * 0.1 * 0.0045);                      % Volumetric dissipation [W/m^3]

        %%%%%%%%%%%% CALCULATIONS %%%%%%%%%%%%
        %%%% SECTION 1: From the PCB border to the start of the 1st IC 
        T = zeros(1, M);
        for i = 1: lim(1)
            T(i) = temp_lin(Q1, x(i), k_eff, A, a);
        end
        %%%% SECTION 2: From the start of the 1st IC to the end of that IC
        switch(method)
            case 1 %%% k -->  k_ic
                a = T(lim(1));                              % [K]
                b = b_coef(Q1, k_eff_ic, A_ic);             % [K*m]
                for i = lim(1):lim(2)
                    T(i) = temp_parb(a, b, (x(i) - x(lim(1))), phi_ic, k_eff_ic);   % [K]
                end

            case 2 %%% k --> inf
                for i = lim(1):lim(2)
                    T(i) = T(lim(1));                       % [K]
                end
        end
        %%%% SECTION 3: From the end of the 1st IC to the start of the 2nd IC
        for i = (lim(2)+1):lim(3)
            T(i) = temp_lin(Q2, (x(i) - x(lim(2))), k_eff, A, T(lim(2)));           % [K]
        end
        %%%% SECTION 4: From the start of the 2nd IC to the center of that IC
        switch(method)
            case 1 %%% k -->  k_ic
                a = T(lim(3));                              % [K]
                b = b_coef(Q2, k_eff_ic, A_ic);             % [K*m]
                for i = (lim(3)+1): (M/2)
                    T(i) = temp_parb(a, b, (x(i) - x(lim(3))), phi_ic, k_eff_ic);   % [K]
                end

            case 2 %%% k --> inf
                for i = (lim(3)+1): (M/2)
                    T(i) = T(lim(3));                       % [K]
                end
        end
        %%%% TAKE ADVANTAGE OF SYMMETRY
                T(((M/2)+1):M) = T(M/2:-1:1);                       % Mirror curve

                T_0 = max(T)                                        % Max T [K]
                T_0_C = convtemp(T_0, 'K', 'C')                     % Max T [Celsius]
               
        %%%% PLOT TEMPERATURE PROFILE
                figure()
                hold on
                myplot(x,T)
                hold on
                axis([0 dx T_b*0.9 max(T)*1.1])
                ylabel('{\it T} [K]')
                xlabel('{\it x} [m]');
                xline(L, '-.')
                yline(T_b, '--')
                yline(T_0, '--')
                hold off
%___________________________________________________________________________
        %% Apartado C

    case 'c'

%___________________________________________________________________________
        %% Apartado D

    case 'd'

%___________________________________________________________________________
        %% Apartado E

    case 'e'
        
%___________________________________________________________________________
end

%___________________________________________________________________________
%% ======= FUNCIONES ADICIONALES ======= %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective thermal conductivity k
function k_eff = effective(k_vect, l_vect)
k_eff = sum(k_vect.*l_vect)/sum(l_vect);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parabolic temperature eq
function T = temp_parb(a, b, x, phi, k)
T = a + b * (x) - ( phi/(2*k) ) * (x)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'b' parameter in parabolic temperature eq
function b = b_coef(Q, k, A)
b = Q / (k * A);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear termperature eq
function T = temp_lin(Q, x, k, A, a)
DT = Q * x / (k * A);
T = a + DT;
end
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