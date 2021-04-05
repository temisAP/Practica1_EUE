%% Datos
DT = 50; %K
L = 0.25; %m


%% Elementos

% Material
mat = material();
mat.A = 0.0001; %m^2
mat.E = 70e9; %Pa
mat.alpha = 22.5e-6; %K^-1

% Barras
b12 = barra(mat);
b12.L = L/2^0.5; %m
b12.theta = 0;
b12.nodes = [1 2];

b13 = barra(mat);
b13.L = L; %m
b13.theta = deg2rad(45);
b13.nodes = [1 3];

b23 = barra(mat);
b23.L = L/2^0.5; %m
b23.theta = deg2rad(90);
b23.nodes = [2 3];

%% Estructura

bs = [b12,b13,b23];
s = estructura(bs,2);

%% Desplazamientos

% Nodos libres
a = [1,6];
Ka = s.K(a,a);
sigma_t = b13.E * b13.alpha * DT;
Fa = [-cos(b13.theta) sin(b13.theta)]' * sigma_t * b13.A ;

Ua = inv(Ka)*Fa;

%% Reacciones en los apoyos

% Movimientos restringidos
c = [2,3,4,5];

Kca = s.K(c,a);
Fc = Kca * Ua;

F = zeros(length(s.K),1);
F(a) = Fa;
F(c) = Fc;

%% Tensión sobre cada varilla

U=zeros(length(s.K),1);
U(a) = Ua;
U(c) = zeros(size(c));

% Varilla 1-2

u1 = U(1);
u2 = U(3);

b12.F = b12.K * [u1 0 u2 0]';
b12.sigma = (b12.F(3)-b12.F(1))/2 / b12.A;

% Varilla 1-3

u1 = U(1)*cos(b13.theta) + U(2)*sin(b13.theta);
u2 = U(5)*cos(b13.theta) + U(6)*sin(b13.theta);

b13.F = b13.K * [u1 0 u2 0]';
b13.sigma = (b13.F(3)-b13.F(1))/2/ b13.A;

% Varilla 2-3

u1 = U(4);
u2 = U(6);

b23.F = b23.K * [u1 0 u2 0]';
b23.sigma = (b23.F(3)-b23.F(1))/2/ b23.A;

%% Resultados

disp('Fuerzas de reacción en los apoyos')
disp(Fc)

disp('Desplazamiento de los nodos 1 y 3')
disp(Ua)

disp('Tensión sobre cada varilla')
disp(' b12')
disp(b12.sigma)
disp(' b13')
disp(b13.sigma)
disp(' b23')
disp(b23.sigma)

%% Deformada de la estructura

factor = 100;

% Analítica
U_ini = [ -b12.L,0 , 0,0 , 0,b23.L];
U_fin = U_ini + U'*factor;

plot_estructura(U_ini,U_fin,'Figuras/deformada_analitica')

%Numérica
load('numerico/estatico6num.mat')
U_ini = [ -b12.L,0 , 0,0 , 0,b23.L];
U_fin_num = U_ini + U_num'*factor;

plot_estructura(U_ini,U_fin_num,'Figuras/deformada_numerica')

%Las dos
plot_estructura_conjunta(U_ini,U_fin,U_fin_num,'Figuras/deformadasss')

function plot_estructura(U_ini,U_fin,name)

x_ini = [U_ini(1:2:end), U_ini(1)];
y_ini = [U_ini(2:2:end), U_ini(2)];

x_fin = [U_fin(1:2:end), U_fin(1)];
y_fin = [U_fin(2:2:end), U_fin(2)];

    fsz = 10; %Tamaño de la fuente

h=figure();
    hold on
    for i=1:length(x_ini)-1
        plot(x_ini([i,i+1]), y_ini([i,i+1]),'-.', 'LineWidth', 1.5, 'Color','k');
        plot(x_fin([i,i+1]), y_fin([i,i+1]),'-', 'LineWidth', 1.5, 'Color','k');
    end
    plot(x_ini, y_ini,'o', 'LineWidth', 1.5, 'Color','k');
    plot(x_fin, y_fin,'o', 'LineWidth', 1.5, 'Color','k');
    ylabel({'y[m]'}),  xlabel('x[m]'), grid on, box on;

    legend([{'Inicial','Final'}],'Interpreter', 'Latex')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fsz);
    set(gca,'LabelFontSizeMultiplier',1.2);
    set(gca,'TitleFontSizeMultiplier',1.35);
    ylim([-y_fin(3)*0.2 y_fin(3)*1.2]);
    xlim([x_fin(1)*1.2 -x_fin(1)*0.5]);

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
    ylh.Position(1) = -(0.8 + y_lab/100)*( ax(2) - ax(1) );
    ylh.Position(2) = 0.8*( ax(4) - ax(3) ) + ax(3);

    y_pdf = 0;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+y_pdf/10, pos(4)])
    print(h, name,'-dpdf','-r0','-painters')
end

function plot_estructura_conjunta(U_ini,U_fin,U_fin_num,name)

x_ini = [U_ini(1:2:end), U_ini(1)];
y_ini = [U_ini(2:2:end), U_ini(2)];

x_fin = [U_fin(1:2:end), U_fin(1)];
y_fin = [U_fin(2:2:end), U_fin(2)];

x_fin_num = [U_fin_num(1:2:end), U_fin_num(1)];
y_fin_num = [U_fin_num(2:2:end), U_fin_num(2)];

    fsz = 10; %Tamaño de la fuente

h=figure();
    hold on
    for i=1:length(x_ini)-1
        plot(x_ini([i,i+1]), y_ini([i,i+1]),'-', 'LineWidth', 1.5, 'Color','k');
        plot(x_fin([i,i+1]), y_fin([i,i+1]),'-.', 'LineWidth', 1.5, 'Color','c');
        plot(x_fin_num([i,i+1]), y_fin_num([i,i+1]),'--', 'LineWidth', 1.5, 'Color','r');
    end
    plot(x_ini, y_ini,'o', 'LineWidth', 1.5, 'Color','k');
    plot(x_fin, y_fin,'o', 'LineWidth', 1.5, 'Color','k');
    plot(x_fin_num, y_fin_num,'o', 'LineWidth', 1.5, 'Color','k');
    ylabel({'y[m]'}),  xlabel('x[m]'), grid on, box on;

    legend([{"Inicial","Anal\'itica","Num\'erica"}],'Interpreter', 'Latex')
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize',fsz);
    set(gca,'LabelFontSizeMultiplier',1.2);
    set(gca,'TitleFontSizeMultiplier',1.35);
    ylim([-y_fin(3)*0.2 y_fin(3)*1.2]);
    xlim([x_fin(1)*1.2 -x_fin(1)*0.5]);

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
    ylh.Position(1) = -(0.8 + y_lab/100)*( ax(2) - ax(1) );
    ylh.Position(2) = 0.8*( ax(4) - ax(3) ) + ax(3);

    y_pdf = 0;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+y_pdf/10, pos(4)])
    print(h, name,'-dpdf','-r0','-painters')
end
