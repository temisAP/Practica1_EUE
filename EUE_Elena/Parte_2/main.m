clear all
close all
fig = 1;
fig_flag = 1;

%% Datos

Lx = 1;     % m
Ly = 1.25;  % m
H  = 50e-3; % m

Panel = struct();
Panel.A         = Lx*Ly;    % m^2
Panel.thickness = 5e-3;     % m
Panel.rho       = 2700;     % kg/m^3
Panel.E         = 70e9;     % Pa
Panel.nu        = 0.3;
Panel.D         = Panel.E*Panel.thickness^3/12/(1-Panel.nu^2);

Aire = struct();
Aire.rho        = 1.225;    % kg/m^3
Aire.c          = 343;      % m/s
Aire.V          = Panel.A*H;                      
Aire.A          = 2 * Panel.A + 2*(Lx+Ly) * H ;
Aire.L          = 4*(Lx+Ly+H);

%% /////// APARTADO A /////// %%
%% Número de modos por banda

% Frecuencias en bandas de tercio de octava con la siguiente estrcutura:
%   Límite inferior, Frecuencia central, Límite superior, Incremento
try 
    load('Frecuencias.mat','Frecuencias');
catch
    disp('Error al cargar Frecuencias.mat')
end

% Paneles
Panel.modos(:,:) = Frecuencias(:,:);
Panel.modos(:,5) = Panel.A/4/pi * (Panel.rho*Panel.thickness/Panel.D)^0.5;
Panel.modos(:,6) = Panel.modos(:,5) .* Panel.modos(:,4) ;

% Capas de aire
Aire.modos(:,:) = Frecuencias(:,:);
Aire.modos(:,5) = Aire.V/pi/Aire.c * (Frecuencias(:,2).*2*pi ./ Aire.c).^2 +...
    Aire.A/4/Aire.c * (Frecuencias(:,2).*2*pi ./ Aire.c) + ...
    Aire.L/8/Aire.c;
Aire.modos(:,6) = Aire.modos(:,5) .* Aire.modos(:,4);

% Requisito para el SEA

idx = min(find(Panel.modos(:,6)>5));
Panel.f_SEA = Panel.modos(idx,2);

idx = min(find(Aire.modos(:,6)>5));
Aire.f_SEA = Aire.modos(idx,2);

clear idx

SEA_frec = max([Panel.f_SEA, Aire.f_SEA]);

if fig_flag == 1
    figure(fig)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    hold on
    plot(Panel.modos(:,2),Panel.modos(:,6), 'DisplayName', 'Panel')
    plot(Aire.modos(:,2), Aire.modos(:,6),  'DisplayName', 'Aire' )
    yline(5,'--','Color','k','DisplayName', 'N=5')
    xline(Panel.f_SEA,'Color','g', 'DisplayName', 'Panel f SEA')
    xline(Aire.f_SEA,'Color','y', 'DisplayName', 'Aire f SEA')
    grid on
    ylabel('N (Número de modos)')
    xlabel('Frecuencia [Hz]')
    legend('Location','best')
    fig = fig+1;
end


%% /////// APARTADO B /////// %%
%% Factor de pérdidas de acoplamiento

% Frecuencias de la placa
Panel.m   = Panel.rho * Panel.thickness;              % Masa por unidad de superficie
Panel.fc  = (Aire.c^2/2/pi) * (Panel.m/Panel.D)^0.5;  % Frecuencia crítica de la placa
Panel.f11 = (Aire.c^2/4/Panel.fc) * (1/Lx^2+1/Ly^2);  % Primer modo superficial de la placa

% Sigma

A   = Panel.A;
c0  = Aire.c;
P   = 2*(Lx+Ly);
fc  = Panel.fc;
    
for f = 1:length(Frecuencias(:,2))
    frec   = Frecuencias(f,2);
    lambda = (frec/fc)^0.5 ;
    
    if Panel.f11 <= Panel.fc/2
        if frec<Panel.f11
            sigma(f) = 4*A^2*(frec/c0)^2;
        elseif Panel.f11< frec && frec < Panel.fc/2
            A1 = c0*P/A/fc;
            A2 = ((1-lambda^2)*log((1+lambda)/(1-lambda))+2*lambda) / (4*pi^2*(1-lambda^2)^1.5);
            A3 = 2*(2*c0/fc/pi^2)^2;
            A4 = (1-2*lambda^2)/(A*lambda*(1-lambda^2)^0.5);
            sigma(f) = A1 * A2 + A3 * A4;
        elseif Panel.fc/2 < frec && frec < Panel.fc
            A1 = c0*P/A/fc;
            A2 = ((1-lambda^2)*log((1+lambda)/(1-lambda))+2*lambda) / (4*pi^2*(1-lambda^2)^1.5);
            sigma(f) = A1 * A2;
        elseif frec >= Panel.fc
            sigma(f) = 1/(1-fc/frec)^0.5;
        end
    elseif Panel.f11 > Panel.fc/2
        if frec <= Panel.fc
            sigma(f) = 4*A^2*(frec/c0)^2;
        elseif frec > Panel.fc
            sigma(f) = 1/(1-fc/frec)^0.5;
        end
    end
end

clear A1 A2 A3 A4
clear A c0 P fc
clear deltaf deltaW lambda

sigma = sigma';

% Pérdidas panel-aire
Panel.eta = Aire.A*Aire.rho*Aire.c.*sigma./(Panel.rho*Panel.A*Panel.thickness*2*pi.*Frecuencias(:,2));

% Pérdidas aire-panel
Aire.eta  = Panel.eta .* Panel.modos(:,5) ./ Aire.modos(:,5);

% Representación gráfica

if fig_flag == 1
    figure(fig)
    set(gca, 'XScale', 'log')
    hold on
    plot(Frecuencias(:,2), Panel.eta(:), 'DisplayName', 'Panel')
    plot(Frecuencias(:,2), Aire.eta(:),  'DisplayName', 'Aire' )
    grid on
    legend('Location','best')
    ylabel('eta')
    xlabel('Frecuencia [Hz]')
    fig = fig+1;
end

%% /////// APARTADO C /////// %%

% {P_ext(w)} = w * [eta(w)] * {E(w)}

%% /////// APARTADO D /////// %%

%% Ensamblaje (haciendo uso de la clase conjunto de la parte anterior)

% Panel-Aire

eta_pd = 0.015;

Panel.M(1,1,:) = +Panel.eta(:)+eta_pd;
Panel.M(1,2,:) = -Panel.eta(:);
Panel.M(2,1,:) = -Aire.eta(:);
Panel.M(2,2,:) = +Aire.eta(:);

% Aire-Panel

eta_ad = 0.01;

Aire.M(1,1,:) = +Aire.eta(:)+eta_ad;
Aire.M(1,2,:) = -Aire.eta(:);
Aire.M(2,1,:) = -Panel.eta(:);
Aire.M(2,2,:) = +Panel.eta(:);

% Conjuntos
Paneles = [Panel Panel Panel];      
Aires   = [Aire Aire];

% Nodos
Paneles(1).nodes    = [1 2];
Aires(1).nodes      = [2 3];
Paneles(2).nodes    = [3 4];
Aires(2).nodes      = [4 5];

% Ensamblaje
estructura = conjunto({Paneles(1), Paneles(2), Aires(1), Aires(2)},1);  %Hay que quitar el último panel porque no tiene aire detrás
estructura.M(5,5,:) = estructura.M(5,5,:) + eta_pd;                     %Faltaba la disipación del último panel

%% Potencias externas aplicadas

% Valores de la tabla
Frecuencias_list = [16 1e3 1250 1600 2000 2500 3150 4000 5000 1e4];

Potencias_list = [  10  10   20   35   50   80  100  150  100 100;...
                    4.35 4.35 8.7 15.2 21.74 36.96 39.13 45.65 43.47 43.47;...
                    10  10   20   35   50   80  100  150  100 100] ; % W            
                
% Columna con frecuencias (Hz)
for p=1:length(Paneles)
    Paneles(p).p_ext(:,1) = Frecuencias(:,2);
end

% Columna con potencias (W)

f = 1;                
                            
Paneles(1).p_ext(f,2) = 10;
Paneles(2).p_ext(f,2) = 4.35;
Paneles(3).p_ext(f,2) = 10;
                      
for f=2:length(Frecuencias)
    
    if Frecuencias(f,2) > Frecuencias_list(1) && Frecuencias(f,2) <= Frecuencias_list(2) l = 2; end
    if Frecuencias(f,2) > Frecuencias_list(2) && Frecuencias(f,2) <= Frecuencias_list(3) l = 3; end
    if Frecuencias(f,2) > Frecuencias_list(3) && Frecuencias(f,2) <= Frecuencias_list(4) l = 4; end
    if Frecuencias(f,2) > Frecuencias_list(4) && Frecuencias(f,2) <= Frecuencias_list(5) l = 5; end
    if Frecuencias(f,2) > Frecuencias_list(5) && Frecuencias(f,2) <= Frecuencias_list(6) l = 6; end
    if Frecuencias(f,2) > Frecuencias_list(6) && Frecuencias(f,2) <= Frecuencias_list(7) l = 7; end
    if Frecuencias(f,2) > Frecuencias_list(7) && Frecuencias(f,2) <= Frecuencias_list(8) l = 8; end
    if Frecuencias(f,2) > Frecuencias_list(8) && Frecuencias(f,2) <= Frecuencias_list(9) l = 9; end
    if Frecuencias(f,2) > Frecuencias_list(9) && Frecuencias(f,2) <= Frecuencias_list(10) l = 10; end
    
    for p=1:length(Paneles)
        deltaW = (Potencias_list(p,l) - Potencias_list(p,l-1))  / (Frecuencias_list(l) - Frecuencias_list(l-1)) ;
        deltaf = Frecuencias(f,2) - Frecuencias(f-1);
        Paneles(p).p_ext(f,2) = Paneles(p).p_ext(f-1,2) + deltaW*deltaf;
    end

end

clear l Frecuencias_list Potencias_list

% Comprobación visual

if fig_flag == 1
    figure(fig)
    hold on
    for p=1:length(Paneles)
        plot(Paneles(p).p_ext(:,1),Paneles(p).p_ext(:,2), 'DisplayName', ['Panel',num2str(p)])
    end
    grid on
    legend()
    ylabel('P ext [W]')
    xlabel('Frecuencia [Hz]')
    fig = fig+1;
end

% En un único vector

Z = zeros(length(Frecuencias(:,2)),1);
P_ext = [Paneles(1).p_ext(:,2) Z Paneles(2).p_ext(:,2) Z Paneles(3).p_ext(:,2)]';


%% Energías 

idx = find(Frecuencias(:,2)>SEA_frec);

for f = idx(1):idx(end)
    omega = Frecuencias(f,2) * 2 * pi;
    M = estructura.M(:,:,f); % Porque el último panel no tiene aire luego
    P = P_ext(:,f);
    E(:,f) = 1/omega .* inv(M) * P;
end

clear omega M P

% Representación gráfica

if fig_flag == 1
    figure(fig)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    hold on
    display_names = {'Panel 1', 'Aire 1', 'Panel 2', 'Aire 2', 'Panel 3'};
    for p = 1:length(E(:,1))
        FF = Frecuencias(idx,2);
        EE = E(p,idx);
        plot(FF,EE,'DisplayName',display_names{p});
    end
    grid on
    legend()
    ylabel('E [J]')
    xlabel('Frecuencia [Hz]')
    fig = fig+1;
end

clear EE FF


%% /////// APARTADO E /////// %%

%% Velocidad media en los paneles
% E = 1/2 * m * v^2

v(1:3,:) = (2*E([1,3,5],:)./Panel.A/Panel.rho/Panel.thickness).^0.5;


if fig_flag == 1
    figure(fig)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    hold on
    display_names = {'Panel 1', 'Panel 2' 'Panel 3'};
    for p = 1:length(v(:,1))
        FF = Frecuencias(idx,2);
        vv = v(p,idx);
        plot(FF,vv,'DisplayName',display_names{p});
    end
    grid on
    legend()
    ylabel('v [m/s]')
    xlabel('Frecuencia [Hz]')
    fig = fig+1;
end

clear FF vv

%% Presión media en las capas de aire
% E = V * P^2 / (rho*c0^2)

P(1:2,:) = (E([2,4],:)./Aire.V * Aire.rho * Aire.c^2).^0.5;

if fig_flag == 1
    figure(fig)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    hold on
    display_names = {'Aire 1', 'Aire 2'};
    for p = 1:length(P(:,1))
        FF = Frecuencias(idx,2);
        PP = P(p,idx);
        plot(FF,PP,'DisplayName',display_names{p});
    end
    grid on
    legend()
    ylabel('P [Pa]')
    xlabel('Frecuencia [Hz]')
    fig = fig+1;
end

clear FF PP
