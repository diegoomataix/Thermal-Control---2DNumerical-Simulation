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
dx = 140e-3;    % m
dy = 100e-3;    % m
dz = 1.5e-3;    % m
Vol = dx*dy*dz; % m^2
    % Recubrimiento de cobre del PCB por cada lado
t_rec = 50e-6;  % m     % en una de las caras es continuo, y en la otra 
                        % ocupa solo el 10% de la superficie en la cual van 
                        % montados tres IC
    % los lados cortos de la PCB tienen contacto termico con paredes a 25C
T_paredes =  convtemp(25, 'C', 'K'); 
    % los otros dos bordes estan termicamente aislados.
    % para el FR4:
k_plano = 0.5;              % W / ( m * K )
k_traves = k_plano / 2;     % W / ( m * K )
                        
%%%% IC %%%%
dx_ic = 40e-3;  % m
dy_ic = 20e-3;  % m
dz_ic = 3e-3;   % m
Vol_ic = dx_ic * dy_ic * dz_ic; % m^2
    % dispando:
disip_ic = 5;   % W
    % conductividad termica:
k_ic = 50;      % W / ( m * K )
    % capacidad termica:
C_ic = 20;      % J / K
    % separación entre ICs
dist_ic = 20e-3;% m 

%%%% Otros datos %%%%
% calor transmitido por radiacion:
emis_comp = 0.7; % emisividad media por el lado de los componentes
emis_cara = 0.5; % emisividad media por el lado de la cara opuesta
% caja electronica que se puede considerar negra
T_caja = convtemp(45, 'C', 'K'); % T caja electronica (para apartado c)

