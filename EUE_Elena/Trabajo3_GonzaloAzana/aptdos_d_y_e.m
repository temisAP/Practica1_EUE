clc
close all
clear all

% Cogemos los datos del excel
frec = xlsread('Practica3.xlsx',2,'J17:J32');
etas = xlsread('Practica3.xlsx',2,'U17:V32');

Ps = xlsread('Practica3.xlsx',2,'P17:Q32');


%% Aptdo. d - Calculo de la energia

eta_pa = etas(:,1);
eta_ap = etas(:,2);

eta_pd = 1.5e-2;
eta_ad = 1e-2;

P_1 = Ps(:,1);
P_2 = Ps(:,2);

omega = 2*pi*frec;

Etas = zeros(5,5,length(omega));

P = zeros(5,length(omega));
% eta_p3a = eta_p1_a
% eta_a_p3 = eta_a_p1
% P_3 = P_1
for j = 1:length(omega)
    % Subs. 1
    Etas(1,1,j) = eta_pd + eta_pa(j);
    Etas(1,2,j) = - eta_ap(j);
    %Subs. 2
    Etas(2,1,j) = - eta_pa(j);
    Etas(2,2,j) = eta_ad + eta_ap(j)+ eta_ap(j);
    Etas(2,3,j) = - eta_pa(j);
    % Subs. 3
    Etas(3,2,j) = - eta_ap(j);
    Etas(3,3,j) = eta_pd + 2*eta_pa(j);
    Etas(3,4,j) = - eta_ap(j);
    % Subs. 4
    Etas(4,3,j) = - eta_pa(j);
    Etas(4,4,j) = eta_ad + eta_ap(j) + eta_ap(j);
    Etas(4,5,j) = - eta_pa(j);
    % Subs. 5
    Etas(5,4,j) = - eta_ap(j);
    Etas(5,5,j) = eta_pd + eta_pa(j);
    
    P(1,j) = P_1(j);
    P(5,j) = P_1(j);
    P(3,j) = P_2(j);
end

E = zeros(5,length(omega));
for j = 1:length(omega)
%     size(P(j))
%     size(Etas(:,:,j))
   E(:,j) = (omega(j)*Etas(:,:,j))\ P(:,j);
end

E = E';

figure()
loglog(omega/2/pi,[ E(:,1) E(:,2) E(:,3)])
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$E$ [J]','interpreter','latex')
legend('Placas 1 y 3','Aire','Placa 2','Interpreter','latex')

E = E';

v_m = sqrt(E([1 3],:)/(2700*1*1.25*0.005))';
Pres = sqrt(E(2,:)/1.25/1/0.05*1.225*343^2)';

figure()
loglog(frec,v_m)
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$v$ [m/s]','interpreter','latex')
legend('Placas 1 y 3','Placa 2','Interpreter','latex')

figure()
loglog(frec,Pres)
grid on
xlabel('$f$ [Hz]','interpreter','latex')
ylabel('$P$ [Pa]','interpreter','latex')