clear all
syms eta p sigma emiss T_avg k_eff A phi lambda dx c2 c1 x T_box  T_b  
latex(eta == 4*p * sigma * emiss * T_avg^3)         % Auxiliary function for simplifying the ODE
latex(lambda == sqrt( eta / ( (k_eff * A) ) ))                                  % Eigenvalues of the ODE
latex(c2 == (T_b - T_box - ( (phi * A) / eta ) / (1+exp(-lambda *dx)) ))        % Coef. of the ODE
latex(c1 == c2*exp(-lambda*dx))                                                 % Coef. of the ODE
% 
%         % Temperature profile, the expression is found by solving the ODE
%         % and applying the BC.
%% 
syms T
latex(T == c1* exp(lambda * x ) + c2*exp(-lambda * x ) + T_box + ( (phi*A)/eta))  % [K]
%% 
syms k_vect l_vect a b Q DT k
latex(k_eff == sum(k_vect.*l_vect)/sum(l_vect))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parabolic temperature eq
latex(T == a + b * (x) - ( phi/(2*k) ) * (x)^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'b' parameter in parabolic temperature eq
latex(b == Q / (k * A))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear termperature eq

latex(DT == Q * x / (k * A))
latex(T == a + DT)
