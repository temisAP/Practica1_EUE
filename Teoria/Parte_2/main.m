clear all
close all
fig = 1;
fig_flag = 1;
save_flag = 1;

colors = [0, 0.4470, 0.7410;
          [255, 102, 26]/255;
          0, 0.6, 0.7410;
          [255, 200, 26]/255;
          [255,140,0]/255;
          [139,0,139]/255;
          [50,205,50]/255];

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
    CreateFrec
end

% Paneles
Panel.modos(:,:) = Frecuencias(:,:);
Panel.modos(:,5) = Panel.A/4/pi * (Panel.rho*Panel.thickness/Panel.D)^0.5; % modos/Hz
Panel.modos(:,6) = Panel.modos(:,5) .* Panel.modos(:,4);

% Capas de aire
Aire.modos(:,:) = Frecuencias(:,:);
Aire.modos(:,5) = Aire.V/pi/Aire.c * (Frecuencias(:,2).*2*pi ./ Aire.c).^2 +...
    Aire.A/4/Aire.c * (Frecuencias(:,2).*2*pi ./ Aire.c) + ...
    Aire.L/8/Aire.c; % modos/Hz
Aire.modos(:,6) = Aire.modos(:,5) .* Aire.modos(:,4);

% Requisito para el SEA

idx = min(find(Panel.modos(:,6)>5));
Panel.f_SEA = Panel.modos(idx,2);

idx = min(find(Aire.modos(:,6)>5));
Aire.f_SEA = Aire.modos(idx,2);

clear idx

SEA_frec = max([Panel.f_SEA, Aire.f_SEA]);

if fig_flag == 1
    h = figure(fig); %set(h, 'Visible', 'off')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(get(gca,'ylabel'),'rotation',0);
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',10.5);
    set(gca,'TitleFontSizeMultiplier',1.25);
    set(gca,'LabelFontSizeMultiplier',1.3);
    ylh = get(gca,'ylabel'); ylh.Position(1) = 0.35; ylh.Position(2) = 100;
    xlh = get(gca,'xlabel'); xlh.Position(1) = 2000; xlh.Position(2) = 0.004;
    
    hold on
        plot(Panel.modos(:,2),Panel.modos(:,6),'Color', colors(1,:),'DisplayName', 'Panel')
        plot(Aire.modos(:,2), Aire.modos(:,6), 'Color', colors(2,:), 'DisplayName', 'Aire' )
        yline(5,'--','Color','k','DisplayName', 'N=5')
        xline(Panel.f_SEA,'--','Color', colors(1,:), 'DisplayName', 'Panel f SEA')
        xline(Aire.f_SEA,'--','Color', colors(2,:), 'DisplayName', 'Aire f SEA')
    grid on; box on;
    legend('Interpreter', 'Latex', 'Location', 'Best')
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$N$'},'Interpreter','latex');
    
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/Modos_por_banda.pdf','-dpdf','-r0','-painters')
    end
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
Panel.eta = Panel.A*Aire.rho*Aire.c.*sigma./(Panel.rho*Panel.A*Panel.thickness*2*pi.*Frecuencias(:,2));
eta_pd = 0.015;

% Pérdidas aire-panel
Aire.eta  = Panel.eta .* Panel.modos(:,5) ./ Aire.modos(:,5);
eta_ad = 0.01;

% Representación gráfica

if fig_flag == 1
    h = figure(fig); %set(h, 'Visible', 'off')
    
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')
    set(get(gca,'ylabel'),'rotation',0);
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',10.5);
    set(gca,'TitleFontSizeMultiplier',1.25);
    set(gca,'LabelFontSizeMultiplier',1.3);
    ylh = get(gca,'ylabel'); ylh.Position(1) = 0.35; ylh.Position(2) = 100;
    xlh = get(gca,'xlabel'); xlh.Position(1) = 2000; xlh.Position(2) = -0.00077;
    
    hold on
        plot(Frecuencias(:,2), Panel.eta(:),'Color', colors(1,:),'DisplayName', '$\eta_{pa}$')
        plot(Frecuencias(:,2), Aire.eta(:), 'Color', colors(2,:), 'DisplayName', '$\eta_{ap}$' )
        yline(eta_pd,'Color',colors(3,:),'DisplayName', '$\eta_{pd}$')
        yline(eta_ad,'Color',colors(4,:),'DisplayName', '$\eta_{ad}$')
    grid on; box on;
    legend('Interpreter', 'Latex', 'Location', 'Best')
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
   
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/Factor_de_perdidas.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

%% /////// APARTADO C /////// %%

% {P_ext(w)} = w * [eta(w)] * {E(w)}

%% /////// APARTADO D /////// %%

%% Ensamblaje (haciendo uso de la clase conjunto de la parte anterior)

% Panel-Aire

Panel.M(1,1,:) = +Panel.eta(:)+eta_pd;
Panel.M(2,1,:) = -Panel.eta(:);
Panel.M(1,2,:) = -Aire.eta(:);
Panel.M(2,2,:) = +Aire.eta(:);

% Aire-Panel

Aire.M(1,1,:) = +Aire.eta(:)+eta_ad;
Aire.M(2,1,:) = -Aire.eta(:);
Aire.M(1,2,:) = -Panel.eta(:);
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
     
    if Frecuencias(f,2) > Frecuencias_list(1) && Frecuencias(f,2) <= Frecuencias_list(2) l = 1; end
    if Frecuencias(f,2) > Frecuencias_list(2) && Frecuencias(f,2) <= Frecuencias_list(3) l = 2; end
    if Frecuencias(f,2) > Frecuencias_list(3) && Frecuencias(f,2) <= Frecuencias_list(4) l = 3; end
    if Frecuencias(f,2) > Frecuencias_list(4) && Frecuencias(f,2) <= Frecuencias_list(5) l = 4; end
    if Frecuencias(f,2) > Frecuencias_list(5) && Frecuencias(f,2) <= Frecuencias_list(6) l = 5; end
    if Frecuencias(f,2) > Frecuencias_list(6) && Frecuencias(f,2) <= Frecuencias_list(7) l = 6; end
    if Frecuencias(f,2) > Frecuencias_list(7) && Frecuencias(f,2) <= Frecuencias_list(8) l = 7; end
    if Frecuencias(f,2) > Frecuencias_list(8) && Frecuencias(f,2) <= Frecuencias_list(9) l = 8; end
    if Frecuencias(f,2) > Frecuencias_list(9) && Frecuencias(f,2) <= Frecuencias_list(10) l = 9; end
    
    for p=1:length(Paneles)
        deltaW = (Potencias_list(p,l+1) - Potencias_list(p,l))  / (Frecuencias_list(l+1) - Frecuencias_list(l)) ;
        deltaf = Frecuencias(f,2) - Frecuencias_list(l);
        Paneles(p).p_ext(f,2) = Paneles(p).p_ext(f-1,2) + deltaW*deltaf;
        dW(l) = deltaW;
        df(l) = deltaf;
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
P_ext = P_ext';

%% Energías 

idx = find(Frecuencias(:,2)>SEA_frec);

for f = idx(1):idx(end)
    omega = Frecuencias(f,2) * 2 * pi;
    M = estructura.M(:,:,f); % Porque el último panel no tiene aire luego
    P = P_ext(f,:)';
    E(f,:) = 1/omega .* inv(M) * P;
end

clear omega M P

% Representación gráfica

if fig_flag == 1
    h = figure(fig); %set(h, 'Visible', 'off')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',10.5);
    set(gca,'TitleFontSizeMultiplier',1.25);
    set(gca,'LabelFontSizeMultiplier',1.3);
    

    hold on
    display_names = {'Panel 1', 'Aire 1', 'Panel 2', 'Aire 2', 'Panel 3'};
    for p = 1:length(E(1,:))
        FF = Frecuencias(idx,2);
        EE = E(idx,p);
        plot(FF,EE,'DisplayName',display_names{p});
    end
    grid on; box on;
    legend('Interpreter', 'Latex', 'Location', 'Best')
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$E$';'[J]'},'Interpreter','latex');
    
    % Get axis size
    ax = axis();
    
    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);
    
    % Move y label
    ylh = get(gca,'ylabel');
    set(get(gca,'ylabel'),'rotation',0)
    ylh.Position(2) = 0.7*( ax(4) - ax(3) ) + ax(3);
      
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/Energia.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

clear EE FF


%% /////// APARTADO E /////// %%

%% Velocidad media en los paneles
% E = 1/2 * m * v^2

vrms(:,1:3) = (2*E(:,[1,3,5])./Panel.A/Panel.rho/Panel.thickness).^0.5;

if fig_flag == 1
    h = figure(fig); %set(h, 'Visible', 'off')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',10.5);
    set(gca,'TitleFontSizeMultiplier',1.25);
    set(gca,'LabelFontSizeMultiplier',1.3);
  
    hold on
    display_names = {'Panel 1', 'Panel 2', 'Panel 3'};
    for p = 1:length(vrms(1,:))
        FF = Frecuencias(idx,2);
        vv = vrms(idx,p);
        plot(FF,vv,'DisplayName',display_names{p});
    end
    grid on; box on;
    legend('Interpreter', 'Latex', 'Location', 'Best')
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$v$';'[m/s]'},'Interpreter','latex');
    
    % Get axis size
    ax = axis();
    
    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);
    
    % Move y label
    ylh = get(gca,'ylabel');
    set(get(gca,'ylabel'),'rotation',0)
    ylh.Position(1) = ylh.Position(1) - 100;
    ylh.Position(2) = 0.7*( ax(4) - ax(3) ) + ax(3);
   
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/v.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

clear FF vv

%% Presión media en las capas de aire
% E = V * P^2 / (rho*c0^2)

Prms(:,1:2) = (E(:,[2,4])./Aire.V * Aire.rho * Aire.c^2).^0.5;

if fig_flag == 1
    h = figure(fig); %set(h, 'Visible', 'off')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',10.5);
    set(gca,'TitleFontSizeMultiplier',1.25);
    set(gca,'LabelFontSizeMultiplier',1.3);
    
    hold on
    display_names = {'Aire 1','Aire 2'};
    for p = 1:length(Prms(1,:))
        FF = Frecuencias(idx,2);
        PP = Prms(idx,p);
        plot(FF,PP,'DisplayName',display_names{p});
    end
    grid on; box on;
    legend('Interpreter', 'Latex', 'Location', 'Best')
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$P$';'[Pa]'},'Interpreter','latex');
    
    % Get axis size
    ax = axis();
    
    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);
    
    % Move y label
    ylh = get(gca,'ylabel');
    set(get(gca,'ylabel'),'rotation',0)
    ylh.Position(1) = ylh.Position(1) - 100;
    ylh.Position(2) = 0.7*( ax(4) - ax(3) ) + ax(3);
    

    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/Prms.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

clear FF PP
