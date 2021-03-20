clear all; clear variables;clc
Conduccion_NumSim_DATOS
%%% Mesh %%%

N = 9e5;                                                    % # of time steps
%N = 2e6;                                                    % # of time steps
%tsim = 1400;                                                % Total simulation time [s]
tsim = 4000;                                                % Total simulation time [s]
m = 1e1;                % Spatial Subdivisions
Mx = 14*m;              % Total n of spatial subdivisions (x-direction)
my = 1e1;               % Spatial Subdivisions (y-direction)
My = 12*my;             % Total n of spatial subdivisions (y-direction)
x = linspace(0, dx, Mx);                                     % Spatial Coordinates [m]
%%% Initialise %%%
Dx=dx/Mx;               % Element width (x-direction)
Dy=dy/My;               % Element width (y-direction)
X=linspace(0,dx,Mx+1);  % Node position list (equispaced) (x-direction)
Y=linspace(0,dy,My+1);  % Node position list (equispaced) (y-direction)
Dt=tsim/N;              % Time step (you might fix it instead of tsim)
t=linspace(0,tsim,N)';  % Time vector

for i = 1:6
    limx(i) = 2*i*m;
end
for i = 1:10
    limy(i) = 1*i*my;
end

%%%%% x-direction %%%%%
%%% NO IC SEGMENTS %%%
for i = 1:(Mx/2)
    for k = 1:limy(4)
        phixy(i,k) = 0;
    end
end
for i = 1:limx(1)
    for k = limy(4):(My/2)
        phixy(i,k) = 0;
    end
end
for i = limx(2)+1:limx(3)
    for k = limy(4):(My/2)
        phixy(i,k) = 0;
    end
end
%%% IC SEGMENTS %%%
for i = limx(1)+1:limx(2)
    for k = limy(4):(My/2)
        phixy(i,k) = 10;
    end
end
for i = limx(3)+1:(Mx/2)
    for k = limy(4):(My/2)
        phixy(i,k) = 10;
    end
end

surf(phixy)

phixy( ( ( (Mx/2)+1):Mx), 1:(My/2)) =   phixy( ( ( (Mx/2)):-1:1), 1:(My/2));
phixy( ( 1:(Mx)),((My/2)+1):My) =   phixy( ( 1:(Mx)),((My/2)):-1:1);

phixy(Mx+1,:) = phixy(Mx,:);
phixy(:,My+1) = phixy(:,My);

% surf(phixy)
surf(Y,X,phixy)