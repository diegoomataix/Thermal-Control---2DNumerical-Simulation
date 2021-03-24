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






%%%%%%%%%%%%%%%%%% AUX %%%%%%%%%%%%%%%%%%
%         %%% Bidimensional temperature profile equation by means of finite elements methods %%%
%         j=1; T(j,:,:)=T_b;       % Initial temperature profile T(x,t)=0 (assumed uniform)
% %         T(:,Mx+1,:)=T_b;
%         for j = 2:N              % 
%             for i = 2:Mx
%                 for k = 2:My
% %                     Dtrcz(i,k)=(Dt*V2d(i,k))/(C(i,k)*z(i,k));
% %                     kzLaplx(j,i,k)=((k_effxy(i,k)+k_effxy(i+1,k))/2)*((z(i,k)+z(i+1,k))/2)*((T(j-1,i+1,k) - ...
% %                                     T(j-1,i,k))/(Dx^2)) - ((k_effxy(i,k)+k_effxy(i-1,k))/2)*((z(i,k)+z(i-1,k))/2)*...
% %                                     ((T(j-1,i,k) - T(j-1,i-1,k))/(Dx^2));
% %                     kzLaply(j,i,k)=((k_effxy(i,k)+k_effxy(i,k+1))/2)...
% %                                     *((z(i,k) + z(i,k+1))/2)*((T(j-1,i,k+1)-T(j-1,i,k))/(Dy^2)) -...
% %                                     ((k_effxy(i,k)+k_effxy(i,k-1))/2)*((z(i,k)+z(i,k-1))/2)*((T(j-1,i,k) -...
% %                                     T(j-1,i,k-1))/(Dy^2));
% %                     hDT(j,i,k)=(emiss_cara+emiss_comp)*stefan_boltz*(T(j-1,i,k)^4 - T_box^4);
% % 
% %                     T(j,i,k) = T(j-1,i,k) + (Dtrcz(i,k))*...
% %                         ((kzLaplx(j,i,k)) + (kzLaply(j,i,k)) + phi2d(i,k)*z(i,k) - (hDT(j,i,k)));
% 
% %                     Dtrcz(i,k)=(Dt*V2d(i,k))/(C(i,k)*z(i,k));
% %                     kzLaplx(j,i,k)=((k_effxy(i,k)+k_effxy(i+1,k))/2)*((z(i,k)+z(i+1,k))/2)*((T(j-1,i+1,k) - ...
% %                                     T(j-1,i,k))/(Dx^2)) - ((k_effxy(i,k)+k_effxy(i-1,k))/2)*((z(i,k)+z(i-1,k))/2)*...
% %                                     ((T(j-1,i,k) - T(j-1,i-1,k))/(Dx^2));
% %                     kzLaply(j,i,k)=((k_effxy(i,k)+k_effxy(i,k+1))/2)...
% %                                     *((z(i,k) + z(i,k+1))/2)*((T(j-1,i,k+1)-T(j-1,i,k))/(Dy^2)) -...
% %                                     ((k_effxy(i,k)+k_effxy(i,k-1))/2)*((z(i,k)+z(i,k-1))/2)*((T(j-1,i,k) -...
% %                                     T(j-1,i,k-1))/(Dy^2));
% %                     hDT(j,i,k)=(emiss_cara+emiss_comp)*stefan_boltz*(T(j-1,i,k)^4 - T_box^4);
% % 
% %                     T(j,i,k) = T(j-1,i,k) + ((Dt*V2d(i,k))/(C(i,k)*z(i,k)))*...
% %                             (((k_effxy(i,k)+k_effxy(i+1,k))/2)*((z(i,k)+z(i+1,k))/2)*((T(j-1,i+1,k) - ...
% %                                     T(j-1,i,k))/(Dx^2)) - ((k_effxy(i,k)+k_effxy(i-1,k))/2)*((z(i,k)+z(i-1,k))/2)*...
% %                                     ((T(j-1,i,k) - T(j-1,i-1,k))/(Dx^2)) + ((k_effxy(i,k)+k_effxy(i,k+1))/2)...
% %                                     *((z(i,k) + z(i,k+1))/2)*((T(j-1,i,k+1)-T(j-1,i,k))/(Dy^2)) -...
% %                                     ((k_effxy(i,k)+k_effxy(i,k-1))/2)*((z(i,k)+z(i,k-1))/2)*((T(j-1,i,k) -...
% %                                     T(j-1,i,k-1))/(Dy^2)) + phi2d(i,k)*z(i,k) - ((emiss_cara+emiss_comp)*...
% %                                     stefan_boltz*(T(j-1,i,k)^4 - T_box^4)));
% %                                 
% %                                 Tdistr5(k,i,j) = Tdistr5(k,i,j-1) + ((deltaT5*Vxy(k,i))/(C_effxy(k,i)*zxy(k,i)))*...
% %                                     (((k_effxy(k,i+1)+k_effxy(k,i))/2)*((zxy(k,i+1)+zxy(k,i))/2)*(Tdistr5(k,i+1,j-1) - ...
% %                                     Tdistr5(k,i,j-1))/(deltaX5^2) - ((k_effxy(k,i)+k_effxy(k,i-1))/2)*((zxy(k,i)+zxy(k,i-1))/2)*...
% %                                     (Tdistr5(k,i,j-1) - Tdistr5(k,i-1,j-1))/(deltaX5^2) + ((k_effxy(k+1,i)+k_effxy(k,i))/2)...
% %                                     *((zxy(k+1,i)+zxy(k,i))/2)*(Tdistr5(k+1,i,j-1) - Tdistr5(k,i,j-1))/(deltaY5^2) - ...
% %                                     ((k_effxy(k,i)+k_effxy(k-1,i))/2)*((zxy(k,i)+zxy(k-1,i))/2)*(Tdistr5(k,i,j-1) - ...
% %                                     Tdistr5(k-1,i,j-1))/(deltaY5^2) + phixy(k,i)*zxy(k,i) - (eps_sup + eps_inf)*...
% %                                     sigma*(Tdistr5(k,i,j-1)^4 - Tinf^4));
% %                                 
%                                 T(j,i,k) = T(j-1,i,k) + ((Dt*V2d(i,k))/(C(i,k)*z(i,k)))*...
%                                     (((k_effxy(i+1,k)+k_effxy(i,k))/2)*((z(i+1,k)+z(i,k))/2)*(T(j-1,i+1,k) - ...
%                                     T(j-1,i,k))/(Dx^2) - ((k_effxy(i,k)+k_effxy(i-1,k))/2)*((z(i,k)+z(i-1,k))/2)*...
%                                     (T(j-1,i,k) - T(j-1,i-1,k))/(Dx^2) + ((k_effxy(i,k+1)+k_effxy(i,k))/2)...
%                                     *((z(i,k+1)+z(i,k))/2)*(T(j-1,i,k+1) - T(j-1,i,k))/(Dy^2) - ...
%                                     ((k_effxy(i,k)+k_effxy(i,k-1))/2)*((z(i,k)+z(i,k-1))/2)*(T(j-1,i,k) - ...
%                                     T(j-1,i,k-1))/(Dy^2) + phi2d(i,k)*z(i,k) - (emiss_cara+emiss_comp)*...
%                                     stefan_boltz*(T(j-1,i,k)^4 - T_box^4));
%                                 
%                                 
% %                     Dtrcz(i,k)=(Dt*V2d(i,k))/(C(i,k)*z(i,k));
% %                     kzLaplx(j,i,k)=((k_effxy(i,k)+k_effxy(i+1,k))/2)*((z(i,k)+z(i+1,k))/2)*((T(j,i+1,k) - ...
% %                                     T(j,i,k))/(Dx^2)) - ((k_effxy(i,k)+k_effxy(i-1,k))/2)*((z(i,k)+z(i-1,k))/2)*...
% %                                     ((T(j,i,k) - T(j,i-1,k))/(Dx^2));
% %                     kzLaply(j,i,k)=((k_effxy(i,k)+k_effxy(i,k+1))/2)...
% %                                     *((z(i,k) + z(i,k+1))/2)*((T(j,i,k+1)-T(j,i,k))/(Dy^2)) -...
% %                                     ((k_effxy(i,k)+k_effxy(i,k-1))/2)*((z(i,k)+z(i,k-1))/2)*((T(j,i,k) -...
% %                                     T(j,i,k-1))/(Dy^2));
% %                     hDT(j,i,k)=(emiss_cara+emiss_comp)*stefan_boltz*(T(j,i,k)^4 - T_box^4);
% % 
% %                     T(j,i,k) = T(j-1,i,k) + (Dtrcz(i,k))*...
% %                         ((kzLaplx(j-1,i,k)) + (kzLaply(j-1,i,k)) + phi2d(i,k)*z(i,k) - (hDT(j-1,i,k)));
% 
% %                     Dtrcz(i,k)= (Dt*V2d(i,k)) / (C(i,k)*z(i,k));
% %                     kzLaplx(j,i,k)=((k_effxy(i,k)+k_effxy(i+1,k))/2)*((z(i,k)+z(i+1,k))/2)*...
% %                                     ( (T(j-1,i+1,k)-T(j-1,i,k) ) / (Dx^2) ) - ...
% %                                     ((k_effxy(i,k)+k_effxy(i-1,k))/2)*((z(i,k)+z(i-1,k))/2)*...
% %                                     ( (T(j-1,i,k)-T(j-1,i-1,k) ) / (Dx^2) );
% %                     kzLaply(j,i,k)=((k_effxy(i,k)+k_effxy(i,k+1))/2)*((z(i,k)+z(i,k+1))/2)*...
% %                                     ( (T(j-1,i,k+1)-T(j-1,i,k) ) / (Dy^2) ) - ...
% %                                     ((k_effxy(i,k)+k_effxy(i,k-1))/2)*((z(i,k)+z(i,k-1))/2)*...
% %                                     ( (T(j-1,i,k)-T(j-1,i,k-1) ) / (Dy^2) );
% %                     hDT(j,i,k)=(emiss_cara+emiss_comp)*stefan_boltz* ((T(j-1,i,k))^4 - (T_box)^4);
% %                     T(j,i,k)= T(j-1,i,k) + Dtrcz(i,k) * ((kzLaplx(j,i,k))+...
% %                                 kzLaply(j,i,k)+phi2d(j-1,i,k)
%                     %Boundory condition in x-nodes 0 and Nx+1:
%                     T(j,1,k) = T_b;
%                     T(j,Mx+1,k) = T_b;
%                     %Boundory condition in y-nodes 0 and Ny+1: (Adiabatic condition)
%                     T(j,i,1) = T(j,i,2);
%                     T(j,i,My+1) = T(j,i,My);
%                 end
%             end
%         end