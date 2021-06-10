%%________________________________________________________________________%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    HIGH FREQUENCY ANALYSIS (Statistical Energy Analysis Method (SEA) )    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Diego Mataix Caballero
%%________________________________________________________________________%

clc; clear all, close all;
%% Datos

% PLACAS
rho_AL = 2700; % kg/m3
L_1 = 1;  % m
L_2 = 1.25; % m
t = 0.005; % m
v_AL= 0.3;
A = L_1 * L_2;
P = 2*L_1 + 2*L_2;
E = 70e9;
M = A*t*rho_AL;
m= M/A;

%AIRE
rho_air = 1.225; % kg/m3
c_0 = 343; %m/sç
h=0.05; % Distancia entre paneles
Lprima = 2*P+4*h;
Aprima = 2*A + 2*(L_1*h) + 2*(L_2*h);
V = A*h;

% Bandas de tercios de octava
f_vector = [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000]; % centro de octava (frecuencia central de la banda de octava)
f_range  = [(17.8-14.1),(22.4-17.8), (28.2-22.4), (35.5-28.2), (44.7-35.5), (56.2-44.7), (70.8-56.2), (89.1-70.8), (112-89.1), (141-112), (178-141), ...
            (224-178), (282-224), (355-282), (447-355), (562-447), (708-562), (891-708), (1122-891), (1413-1122), (1778-1413),...
            (2239-1778),(2818-2239),(3548-2818),(4467-3548),(5623-4467),(7079-5623),(8913-7079),(11220-8913)];

omega = f_vector*2*pi;

% Vectores de potencias externas
Pot_1 = zeros(1,length(f_vector));
Pot_2 = zeros(1,length(f_vector));
Pot_3 = zeros(1,length(f_vector));

Pot_1(1:19) = 10;
Pot_1(20:25)= [20, 35, 50, 80, 100, 150];
Pot_1(26:29) = 100;

Pot_2(1:19) = 4.35;
Pot_2(20:25)= [8.70, 15.20, 21.74, 36.96,39.13, 45.65];
Pot_2(26:29) = 43.47;

Pot_3 = Pot_1;

%% APARTADO 1

% Calculo de la Rigidez de las PLACAS
D = (E*t^3)/(12*(1-v_AL^2));

% Densidad modal del panel
np = (A/(4*pi))*sqrt(rho_AL*t/D);
np_modos = np*f_range;

% Densidad modal del aire
for ii = 1: length(omega)
    na(ii) = ( V / (pi*c_0) ) * ((omega(ii)) / c_0)^2 + (Aprima / (4*c_0)) * ((omega(ii)) / c_0) +  (Lprima / (8*c_0));
end

na_modos = f_range .* na;

% frecuencia crítica de la placa (fc)
fc = (c_0^2/(2*pi))*sqrt(m/D);
lambda = sqrt(f_vector/fc);

% frecuencia del primer modo superficial de la placa (f11)
f11 = (c_0^2/(4*fc)) * ((1/(L_1^2))+(1/(L_2^2)));


% Calculo de la eficiencia de la radiacion
for ii = 1:length(f_vector)
if f11 <= fc/2
    if f_vector(ii) < f11
        sigma(ii) = 4*A^2*(f_vector(ii)/c_0)^2 ;
    elseif f_vector(ii) < (fc/2)
        sigma(ii) = ((c_0*P)/(A*fc)) * ((1-lambda(ii)^2)*log((1+lambda(ii))/(1-lambda(ii))) + 2*lambda(ii)) / ((4*pi^2)*(1-lambda(ii)^2)^(3/2)) ...
                + 2 * ((2*c_0)/(fc*pi^2))^2 * ((1-2*lambda(ii)^2) / (A*lambda(ii)*sqrt(1-lambda(ii)^2)));
    elseif f_vector(ii) < fc
        sigma(ii) =((c_0*P)/(A*fc)) * ((1-lambda(ii)^2)*log((1+lambda(ii))/(1-lambda(ii))) + 2*lambda(ii)) / ((4*pi^2)*(1-lambda(ii)^2)^(3/2));
    else
        sigma(ii) = 1 / (sqrt(1-(fc/f_vector(ii))));
    end
elseif  f11 > fc/2
    if f_vector(ii) < fc
        sigma(ii) = 4*A^2*(f_vector(ii)/c_0)^2 ;
    else
        sigma(ii) = 1 / (sqrt(1-(fc/f_vector(ii))));
    end
end
end

%%
% Esto va a ser un bucle en frecuencias
for ii = 1:length(f_vector)
    n_pa(ii) = (A*rho_air*c_0*sigma(ii))/(M*omega(ii));
end

for ii = 1:length(f_vector)
    n_ap(ii) = (np/na(ii)) * n_pa(ii);
end
% %  % %
% figure()
% plot(f_vector, n_pa)
% hold on
% plot(f_vector, n_ap)
%
% figure()
% loglog(f_vector, n_pa)
% hold on
% loglog(f_vector, n_ap)

%% APARTADO 2
% SISTEMA DE ECUACIONES MATLAB
% subsystem1 = panel; % subsystem2 = aire; % subsystem3 = panel
% subsystem4 = aire; % subsystem5 = panel

n_pd = 0.015;
n_ad = 0.01;

%esta
for i = 22:length(f_range)
Pot_Matrix = [Pot_1(i), 0, Pot_2(i), 0, Pot_3(i)];
n_Matrix = [n_pa(i)+n_pd, -n_ap(i), 0, 0, 0;...
            -n_pa(i), n_ap(i)+n_ad, -n_pa(i), 0, 0;...
            0, -n_ap(i), n_pa(i)+n_pd, -n_ap(i), 0;...
            0, 0, -n_pa(i), n_ap(i)+n_ad, -n_pa(i);...
            0, 0, 0, -n_ap(i), n_pa(i)+n_pd];

B = omega(i).*n_Matrix;
X = linsolve(B,Pot_Matrix');
E_sols (:, i-21) = X;
end

% Velocidad media
M_vector = [M,V*rho_air,M,V*rho_air,M];
v=zeros(size(E_sols));
for ii =1:length(E_sols)
v(:,ii) = sqrt(E_sols(:,ii)./M_vector');
end

for ii =1:length(E_sols)
P_rms(ii) = sqrt(E_sols(2, ii)*rho_air*c_0^2 / V);
end

figure(); hold on;
for i = 1:length(n_Matrix)
plot(f_vector(22:end),v(i,:))
end
grid on;
legend('1','2','3','4','5')
set(gca, 'XScale', 'log')

figure(); hold on;
plot(f_vector(22:end),P_rms)
grid on;
set(gca, 'XScale', 'log')

figure(); hold on;
plot(f_vector(22:end),P_rms)
grid on;
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% for i = 1 : length(f_vector)
% syms E_Matrix [5,1]
% Pot_Matrix = [Pot_1(i), 0, Pot_2(i), 0, Pot_3(i)];
% n_Matrix = [n_pa(i)+n_pd, -n_ap(i), 0, 0, 0;...
%             -n_pa(i), n_ap(i)+n_ad, -n_pa(i), 0, 0;...
%             0, -n_ap(i), n_pa(i)+n_pd, -n_ap(i), 0;...
%             0, 0, -n_pa(i), n_ap(i)+n_ad, -n_pa(i);...
%             0, 0, 0, -n_ap(i), n_pa(i)+n_pd];
% eqns = Pot_Matrix == omega(i).*n_Matrix*E_Matrix;
% S(:,i)= solve(eqns, [E_Matrix])
% En(:,i) = double(S.E_Matrix)
% end
%
% for i = 1 : length(f_vector)
% Pot_Matrix = [Pot_1(i), 0, Pot_2(i), 0, Pot_3(i)];
% syms E1_sym E2_sym E3_sym E4_sym E5_sym
%     eq1 = omega(i)*n_pa(i)*E1_sym - omega(i)*n_ap(i)*E2_sym == Pot_1(i);
%     eq2 = - omega(i)*n_pa(i)*E1_sym + omega(i)*n_ap(i)*E2_sym - n_pa(i)*omega(i)*E3_sym == 0;
%     eq3 = - omega(i)*n_ap(i)*E2_sym + omega(i)*n_pa(i)*E3_sym - n_ap(i)*omega(i)*E4_sym == Pot_2(i);
%     eq4 = - omega(i)*n_pa(i)*E3_sym + omega(i)*n_ap(i)*E4_sym - n_pa(i)*omega(i)*E5_sym == 0;
%     eq5 = - omega(i)*n_ap(i)*E4_sym + omega(i)*n_pa(i)*E5_sym == Pot_3(i);
% eqns = [eq1, eq2, eq3, eq4, eq5];
% assume(E1_sym>=0)
% assume(E2_sym>=0)
% % assume(E3_sym>=0)
% % assume(E4_sym>=0)
% % assume(E5_sym>=0)
% S = solve(eqns,[E1_sym E2_sym E3_sym E4_sym E5_sym]);
% E1(i) = double(S.E1_sym);
% E2(i) = double(S.E2_sym);
% E3(i) = double(S.E3_sym);
% E4(i) = double(S.E4_sym);
% E5(i) = double(S.E5_sym);
% end
% Energy = [E1;E2;E3;E4;E5];
