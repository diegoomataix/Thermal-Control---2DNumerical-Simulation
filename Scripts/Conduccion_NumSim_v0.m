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

%% Datos
Conduccion_NumSim_DATOS

%% Choose exercise to run
choose = 'b';

switch(choose)
    case 'a'
%% Apartado A
% Considerando que la tarjeta sólo evacua calor por los bordes,
% determinar la temperatura máxima que se alcanzaría si toda la disipación
% estuviese uniformemente repartida en la PCB y los IC no influyeran.
        phi = (3 * Q_ic) / Vol                              % Volumetric dissipation [W/m^3]
        e =      [t_rec dz_pcb t_rec];                      % Dimension Vector [m]
        k_vect = [k_Cu k_plano (0.1*k_Cu+0.9*k_plano)];     % Conductivity Vector [W/(m·K)] tercera capa es donde van los IC, cubierta solo al 10% de cobre
        k_eff = sum(k_vect.*e)/sum(e)                       % Effective Conductivity [W/(m·K)]

method = 2;

switch(method)
    case 1
        DT = 1/8 * ( phi * dx^2 / k_eff );                  % Delta T [K]
        T_0 = T_b + DT                                      % Max T [K]
        T_0_C = convtemp(T_0, 'K', 'C')                     % Max T [C]
    case 2
        L = dx/2;                                           % [m]
        b = phi * L / k_eff;                                % [K*m]
        m = 1e6;                                            % Spacial Subdivisions
        M = 14*m;                                           % Total n of spacial subdivisions
        x = linspace(0, dx, M);                             % [m]
        for i = 1:M
            T(i) = T_b + b*x(i) - (phi/(2*k_eff))*x(i)^2;   % [K]
        end
        T_0 = max(T)                                        % Max T [K]
        T_0_C = convtemp(T_0, 'K', 'C')                     % Max T [C]

        figure()
        plot(x,T)
end

%% Apartado B
    case 'b'
% Considerando que la tarjeta sólo evacua calor por los bordes, determinar
% la temperatura máxima que se alcanzaría con un modelo unidimensional en el
% que los IC llegaran hasta los bordes aislados, en el límite kIC→∞, y con la kIC dada.

% Mesh
m = 1e5;                                            % Spatial Subdivisions
M = 14*m;                                           % Total n of spatial subdivisions
x = linspace(0, dx, M);                             % Spatial Coordinates [m]

% Límites de cada tramo
t1 = 2*m;
t2 = 4*m;
t3 = 6*m;
t4 = 8*m;
t5 = 10*m;
t6 = 12*m;

% Define coefficients and some parameters
phi = (3 * Q_ic) / Vol;                             % Volumetric dissipation [W/m^3]
e =      [t_rec dz_pcb t_rec];                      % Dimension Vector [m]
k_vect = [k_Cu k_plano (0.1*k_Cu+0.9*k_plano)];     % Conductivity Vector [W/(m·K)] tercera capa es donde van los IC, cubierta solo al 10% de cobre
k_eff = sum(k_vect.*e)/sum(e);                      % Effective Conductivity [W/(m·K)]
L = dx/2;                                           % [m]
A = 0.1 * 0.0015;                                   % [m^2]
a = T_b;                                            % [K]
b = phi * L / k_eff;                                % [K*m]
Q = 1.5*Q_ic;
Q1 = Q;
Q2 = Q*(1/3);

% e_ic =      [dz_pcb t_rec t_rec 3];                         % Dimension Vector [m]
% k_vect_ic = [k_Cu k_plano (0.1*k_Cu+0.9*k_plano) k_ic];     % Conductivity Vector [W/(m·K)] tercera capa es donde van los IC, cubierta solo al 10% de cobre
% k_eff_ic = ( sum(k_vect_ic.*e_ic)  )/sum(e_ic)              % Effective Conductivity [W/(m·K)]
k_eff_ic = ( k_plano*1.4 + k_Cu*0.05 + (0.1*k_Cu+0.9*k_plano)*0.05+3*k_ic ) /4.5;         % Effective Conductivity [W/(m·K)]
A_ic = 0.0045*0.1;
phi_ic = Q_ic / (0.02 * 0.1 * 0.0045);

% inicializacion
T = zeros(1, 1400);
for i = 1: t1
    DT = Q1 * x(i) / (k_eff * A);
    T(i) = a + DT;
end

a = T(t1);
b = Q1 / (k_eff_ic * A_ic);
for i = t1:t2
    T(i) = a + b * (x(i) - x(t1)) - (phi_ic / (2*k_eff_ic)) * (x(i) - x(t1))^2;
end

for i = (t2+1):t3
    DT = Q2 * (x(i) - x(t2)) / (k_eff*A);
    T(i) = T(t2) + DT;
end

a = T(t3);
b = Q2 / (k_eff_ic*A_ic);
for i = (t3+1):t4
    T(i) = a + b * (x(i) - x(t3)) - (phi_ic / (2*k_eff_ic)) * (x(i) - x(t3) )^2;
end

for i = (t4+1):t5
    DT = -Q2 * (x(i) - x(t4)) / (k_eff * A);
    T(i) = DT + T(t4);
end

a = T(t5);
b = -Q1 / (k_eff_ic*A_ic);
for i = (t5+1):t6
    T(i) = a + b * (x(i) - x(t5)) - (phi_ic / (2*k_eff_ic)) * (x(i) - x(t5))^2;
end

for i = (t6+1):M
    DT = -Q1 * (x(i) - x(t6)) / (k_eff*A);
    T(i) = DT + T(t6);
end

figure()
plot(x,T)
T_0 = max(T)                                        % Max T [K]
T_0_C = convtemp(T_0, 'K', 'C')                     % Max T [C]


%% Apartado C

    case 'c'


%% Apartado D

    case 'd'


%% Apartado E

    case 'e'

end


%% ======= FUNCIONES ADICIONALES ======= %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% COATINGS & SURFACE FINISHES %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = Energy_E(h, v)
E = h * v;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lambda = wavelength(c, v)
lambda = c / v;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total blackbody emissive power
function E_b = tot_blacbody_emiss_p(stefan_boltz, T)
E_b = stefan_boltz * T^4    % total blackbody emissive power
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  spectral emissive power of a blackbody
function E_b = planks_law(C_1, C_2, lambda, T)
E_b = C_1 / ( lambda^5 * exp( C_2 / (lambda*T)) - 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = wiens_displacement(lambda_max_p, T)
T = 2898 / lambda_max_p;    % 2898 micrometer * K
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_lambda_e = emitted_rad_int(delta, Qdot_e, dA, theta, domega, dlambda)
I_lambda_e = delta * Qdot_e / (dA * cos(theta) * domega * dlambda);
% W / ( m^2 * sr * micrometer )
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral hemispherical emissive power
function E_lambda  = spectral_hemispherical_emiss_p(I_lambda_e)
E_lambda = pi * I_lambda_e;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for a blackbody: 0 < epsilon < 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral directional emissivity ( diffuse )
function epsilon_lambda_theta = spectral_dir_emiss(I_lambda_e, I_b_lambda_theta, tipo)
if nargin < 3
    tipo = 'long';
end

switch(tipo) % diapositiva 12
    case 'long'     % calculado con intensidades
        epsilon_lambda_theta = I_lambda_e /  I_b_lambda_theta;
    case 'short'    % igualdad de emissivities
        epsilon_lambda_theta = epsilon_lambda;
    otherwise
        epsilon_lambda_theta = 0;
        disp('Error')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral hemispherical emissivity ( greybody )
function epsilon_lambda = spectral_hemispherical_emiss(E_lambda, E_b_lambda)
epsilon_lambda = E_lambda / E_b_lambda;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total hemispherical emissivity
function epsilon = tot_hemispherical_emiss(E, E_b)
epsilon_lambda = E / E_b;
end
