%% **********************************************************************************
%                      CONDUCCIÓN DE CALOR, SIMULACIÓN NUMÉRICA
% ------------------------------------------------------------------------------------
% Realizado por Diego Mataix Caballero.
%
%  ADDITIONAL NOTES:
%    
%___________________________________________________________________________
close all; clear all; clc;

%% Datos
%%%% PCB (de FR4) %%%%
dx = 140e-3;    % [m]
dy = 100e-3;    % [m]
dz = 1.5e-3;    % [m]
dz_pcb = 1.4e-3;% [m]
dz_cu = (dz - dz_pcb)/2;% [m]
A = dx*dy;      % [m^2]
Vol = dx*dy*dz; % [m^3]
rho_FR4 = 2100; % [kg/m^3]
    % capacidad termica:
c_FR4 = 700;     % [J / kg * K]
c_Cu = 390;      % [J / kg * K]
    % Recubrimiento de cobre del PCB por cada lado
t_rec = 50e-6;  % [m]     % en una de las caras es continuo, y en la otra 
                        % ocupa solo el 10% de la superficie en la cual van 
                        % montados tres IC
rho_Cu = 8960;  % [kg/m^3]
    % los lados cortos de la PCB tienen contacto termico con paredes a 25C
T_b =  convtemp(25, 'C', 'K'); 
    % los otros dos bordes estan termicamente aislados.
    % para el FR4:
k_plano = 0.5;              % [W / ( m * K )]
k_traves = k_plano / 2;     % [W / ( m * K )]
                        
%%%% IC %%%%
dy_ic = 40e-3;  % [m]
dx_ic = 20e-3;  % [m]
dz_ic = 3e-3;   % [m]
Vol_ic = dx_ic * dy_ic * dz_ic; % [m^3]
A_ic = (dz+dz_ic)*dy;           % [m^2]
    % dispando:
Q_ic = 5;                   % [W]
Q_ic_tot = 3 * Q_ic /2;     % [W] (symmetry)
    % conductividad termica:
k_ic = 50;       % [W / ( m * K )]
k_Cu = 395;      % [W / ( m * K )]
    % capacidad termica:
C_ic = 20;       % [J / K]
    % separación entre ICs
dist_ic = 20e-3; % [m] 

%%%% Otros datos %%%%
% calor transmitido por radiacion:
emiss_comp = 0.7; % emisividad media por el lado de los componentes
emiss_cara = 0.5; % emisividad media por el lado de la cara opuesta
% caja electronica que se puede considerar negra
T_box = convtemp(45, 'C', 'K'); % T caja electronica (para apartado c)

%% Constantes

h = 6.6256e-34;  % [J * s]                      % Plank's constant
c_0 = 2.9979e8;  % [m / s]                      % velocidad de la luz
stefan_boltz = 5.67e-8; % [W / ( m^2 * K^4 )]   % Stefan-Boltzmann constant
T_sun = 5800;    % [K]
T_cbr = 2.7;     % [K]
