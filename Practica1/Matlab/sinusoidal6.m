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

P = H*A;                                    %Módulo
theta = atan((2*chi*w/w0) ./ (1-(w/w0).^2));%Argumento

%% Bode de la aceleración

Acc = -w.^2.*P;
Acc = Acc';

h = figure(1);
title('Aceleración')

fsz = 10;

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
ylim([-max(abs(rad2deg(theta)))*1.2 max(abs(rad2deg(theta)))*1.2]);

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
print(h, 'Figuras/Bode','-dpdf','-r0','-painters')

%Save_as_PDF(h,'../Figuras/Bode','horizontal', 7.5, 10);


