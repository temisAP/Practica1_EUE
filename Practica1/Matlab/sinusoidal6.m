%% Datos

% Masa puntual
Mp = 50;     %kg

% Viga
L = 450e-3;     %m
b = 20e-3;      %m (perfil cuadrado)
I = 1/12 * b^4; %m^4
E = 200e9;      %Pa
rho = 7900;     %kg/m^3
nu = 0.3;

% Amortiguamiento
chi = 0.002;

% Excitación
A = 5*9.81; %m/s^2
fmin = 5;   %Hz
fmax = 100; %Hz


%% Otros parámetros

% Rigidez de la viga
% De la deformada bajo una carga unitaria en el centro

Kv = 48*E*I/L^3;

% Masa de la viga

mv = L*b*b * rho;

%% Vibraciones sinusoidales

% Coeficientes (1gdl)
M = Mp + mv;
K = Kv;
F = chi * 2*(K*M)^0.5;

w0 = (K/M)^0.5;
chi = chi;

% Función de transferencia

f = linspace(fmin,fmax,96*2); %Frecuencias discretizadas
w = f*2*pi;

H = (1/K) ./ ((1-(w/w0).^2).^2 + (2*chi*w/w0).^2).^0.5;
theta = atan((2*chi*w/w0) ./ (1-(w/w0).^2));

% Respuesta en frecuencia

P = H*A;                                        %Módulo
theta = atan((2*chi*w/w0) ./ (1-(w/w0).^2));    %Argumento

for i = 1:length(w)
    if w(i)>w0
        theta(i) = theta(i) + deg2rad(180);
    end
end

%% Bode de la aceleración

% Analítico

Acc = -w.^2.*P;
Acc = Acc';

name = 'Figuras/Bode';

el_bode_loco(Acc,theta,w,w0,name)

% Numérico

load('./numerico/sinusoidal6num.mat')

Acc_num = A(:,2);
theta_num = deg2rad(-A(:,3)+180);
w_num = 2*pi*A(:,1);

name = 'Figuras/Bode_num';
el_bode_loco(Acc_num,theta_num,w_num,w0,name)

% Los dos

name = 'Figuras/Bodesss';
el_otro_bode_loco(Acc,theta,w,Acc_num,theta_num,w_num,w0,name)

%% Funciones 

function el_bode_loco(Acc,theta,w,w0,name)

    %Acc es el módulo
    %theta es el argumento

    fsz = 10; %Tamaño de la fuente

    h = figure();
    title('Aceleración')

    subplot(211);
    plot(w/w0, abs(Acc),'-', 'LineWidth', 1.5, 'Color','k');
    title('Módulo'), ylabel({'|Acc|','[m/s^2]'}), xlabel('\Omega/\omega_0'), grid on, box on;
    ylim([-max(abs(Acc))*0.2 max(abs(Acc))*1.2]);

    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fsz);
    set(gca,'LabelFontSizeMultiplier',1.2);
    set(gca,'TitleFontSizeMultiplier',1.35);

    % Get axis size
    ax = axis();

    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);

    % Rotate y label
    set(get(gca,'ylabel'),'rotation',0)
    % Number of characters
    label = get(gca,'ylabel');
    for i = 1:length(label.String)
        ylab_size(i) = length(label.String{i});
    end

    % Move y label
    y_lab = 0;
    ylh = get(gca,'ylabel');
    ylh.Position(1) = -(0.11 + y_lab/100)*( ax(2) - ax(1) );
    ylh.Position(2) = 0.8*( ax(4) - ax(3) ) + ax(3);

    subplot(212);
    plot(w/w0, rad2deg(theta),'-', 'LineWidth', 1.5, 'Color','k');
    title('Argumento'),ylabel({'\theta [º]'}),  xlabel('\Omega/\omega_0'), grid on, box on;
    ini = rad2deg(theta(1));
    fin = rad2deg(theta(end));
    ylim([min([ini,fin])-90 max([ini,fin])+90 ]);

    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fsz);
    set(gca,'LabelFontSizeMultiplier',1.2);
    set(gca,'TitleFontSizeMultiplier',1.35);

    % Get axis size
    ax = axis();

    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);

    % Rotate y label
    set(get(gca,'ylabel'),'rotation',0)
    % Number of characters
    label = get(gca,'ylabel');
    for i = 1:length(label.String)
        ylab_size(i) = length(label.String{i});
    end

    % Move y label
    y_lab = 0;
    ylh = get(gca,'ylabel');
    ylh.Position(1) = -(0.11 + y_lab/100)*( ax(2) - ax(1) );
    ylh.Position(2) = 0.8*( ax(4) - ax(3) ) + ax(3);

    y_pdf = 0;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+y_pdf/10, pos(4)])
    print(h, name,'-dpdf','-r0','-painters')

end

function el_otro_bode_loco(Acc,theta,w,Acc_num,theta_num,w_num,w0,name)

    %Acc es el módulo
    %theta es el argumento

    fsz = 10; %Tamaño de la fuente

    h = figure();
    title('Aceleración')

    subplot(211);
    hold on
    plot(w/w0, abs(Acc),'-.', 'LineWidth', 1.5, 'Color','k','DisplayName',["Anal\'itico"]);
    plot(w_num/w0, abs(Acc_num),'-', 'LineWidth', 1.5, 'Color','k','DisplayName',["Num\'erico"]);
    title('Módulo'), ylabel({'|Acc|','[m/s^2]'}), xlabel('\Omega/\omega_0'), grid on, box on;
    ylim([-max(abs(Acc))*0.2 max(abs(Acc))*1.8]);

    legend('Interpreter', 'Latex')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fsz);
    set(gca,'LabelFontSizeMultiplier',1.2);
    set(gca,'TitleFontSizeMultiplier',1.35);

    % Get axis size
    ax = axis();

    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);

    % Rotate y label
    set(get(gca,'ylabel'),'rotation',0)
    % Number of characters
    label = get(gca,'ylabel');
    for i = 1:length(label.String)
        ylab_size(i) = length(label.String{i});
    end

    % Move y label
    y_lab = 0;
    ylh = get(gca,'ylabel');
    ylh.Position(1) = -(0.11 + y_lab/100)*( ax(2) - ax(1) );
    ylh.Position(2) = 0.8*( ax(4) - ax(3) ) + ax(3);
    hold off

    subplot(212);
    
    hold on
    plot(w/w0, rad2deg(theta),'-.', 'LineWidth', 1.5, 'Color','k','DisplayName',["Anal\'itico"]);
    plot(w_num/w0, rad2deg(theta_num),'-', 'LineWidth', 1.5, 'Color','k','DisplayName',["Num\'erico"]);
    title('Argumento'),ylabel({'\theta [º]'}),  xlabel('\Omega/\omega_0'), grid on, box on;
    ini = rad2deg(theta(1));
    fin = rad2deg(theta(end));
    ylim([min([ini,fin])-90 max([ini,fin])+90 ]);

    legend('Interpreter', 'Latex')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fsz);
    set(gca,'LabelFontSizeMultiplier',1.2);
    set(gca,'TitleFontSizeMultiplier',1.35);

    % Get axis size
    ax = axis();

    % Move x label
    xlh = get(gca,'xlabel');
    xlh.Position(1) = 0.9*( ax(2) - ax(1) ) + ax(1);

    % Rotate y label
    set(get(gca,'ylabel'),'rotation',0)
    % Number of characters
    label = get(gca,'ylabel');
    for i = 1:length(label.String)
        ylab_size(i) = length(label.String{i});
    end

    % Move y label
    y_lab = 0;
    ylh = get(gca,'ylabel');
    ylh.Position(1) = -(0.11 + y_lab/100)*( ax(2) - ax(1) );
    ylh.Position(2) = 0.8*( ax(4) - ax(3) ) + ax(3);

    y_pdf = 0;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+y_pdf/10, pos(4)])
    print(h, name,'-dpdf','-r0','-painters')
    hold off
end



