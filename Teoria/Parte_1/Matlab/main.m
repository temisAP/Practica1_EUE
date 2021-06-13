clear all
close all
fig = 1;
fig_flag = 1;
save_flag = 1;

colors = [0, 0.4470, 0.7410;
          [255, 102, 26]/255;
          [255,140,0]/255;
          [139,0,139]/255;
          [50,205,50]/255];

%% Importar datos de vigas del .pch 

try 
    load('Data.mat')
catch
    ExtractData
end

%% /////// APARTADO A /////// %%
%% Creación de estructura completa

% Renumeración de cada viga
for viga_indx = 1:length(viga)
    viga(viga_indx).N = viga(viga_indx).N - min(viga(viga_indx).N) + 1 ;
end

% Condiciones de contorno
viga(1).bn = [viga(1).N(1)];
viga(1).iface = [viga(1).N(end)];

viga(2).bn = [viga(2).N(end)];
viga(2).iface = [viga(2).N(1)];

% Creación de estructuras 
viga1 = estructura({viga(1)},1);
viga2 = estructura({viga(2)},1);

% Renumeración para que encajen
viga2.N = viga2.N + viga1.N(end) - 1;

% Condiciones de contorno
viga1.bn = [viga1.N(1)];
viga1.iface = [viga1.N(end)];
viga2.bn = [viga2.N(end)];
viga2.iface = [viga2.N(1)];

% Estrcuctura completa
vigas = estructura({viga1,viga2},1);

%% Velocidad rms 

frecuencias = linspace(1,2000,2000); % Hz

for f=1:length(frecuencias)
    q = vigas.dinamic_solution(frecuencias(f));
    qrms(f) = rms(q);
end

vrms = qrms * 2 * pi .*frecuencias;

%% Representación gráfica

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
    xlh = get(gca,'xlabel'); xlh.Position(1) = 2000; xlh.Position(2) = 0.005;
    
    hold on
    plot(frecuencias, vrms, ...
        'LineWidth', 1, 'Color', colors(1,:))
    grid on; box on;
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$v_{rms}$';'[m/s]'},'Interpreter','latex');
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/ApartadoA.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

%% /////// APARTADO B /////// %%
%% Reducción de la estructura

% Condiciones de contorno
viga(1).bn = [viga(1).N(1)];
viga(1).iface = [viga(1).N(end)];

viga(2).bn = [viga(2).N(end)];
viga(2).iface = [viga(2).N(1)];

% Creación de estructuras 
viga1 = estructura({viga(1)},1);
viga2 = estructura({viga(2)},1);

% Transformación Craig-Bampton de las vigas
viga1.CB_transformation();
viga1.CB_getGDLred(5);

viga2.CB_transformation();
viga2.CB_getGDLred(5);

% Renumeración para que encajen
viga2.N = viga2.N + viga1.N(end) - 1;

% Condiciones de contorno
viga1.bn = [viga1.N(1)];
viga1.iface = [viga1.N(end)];
viga2.bn = [viga2.N(end)];
viga2.iface = [viga2.N(1)];

% Estructura completa
vigas_CB = estructura({viga1,viga2},1);

% Calcular las matrices CB a partir de los componentes
vigas_CB.CB_matriz();

% Reducción
vigas_CB.CB_reduction();

%% Velocidad rms

frecuencias = linspace(1,2000,2000); % Hz


for f=1:length(frecuencias)
    q = vigas_CB.CB_sol(frecuencias(f));
    qrms_CB(f) = rms(q);
end

vrms_CB = qrms_CB * 2 * pi .*frecuencias;

%% Representación gráfica

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
    xlh = get(gca,'xlabel'); xlh.Position(1) = 2000; xlh.Position(2) = 0.005;
    
    hold on
    plot(frecuencias, vrms_CB, ...
        'LineWidth', 1, 'Color', colors(1,:))
    
    grid on; box on;
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$v_{rms}$';'[m/s]'},'Interpreter','latex');
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/ApartadoB.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

%% /////// APARTADO C /////// %%


%% Representación gráfica

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
    xlh = get(gca,'xlabel'); xlh.Position(1) = 2000; xlh.Position(2) = 0.005;
    
    hold on
    plot(frecuencias, vrms, ...
        'LineWidth', 1, 'Color', colors(1,:), 'DisplayName', "Soluci\'on completa")
    plot(frecuencias, vrms_CB, ...
        'LineWidth', 1, 'Color', colors(2,:), 'DisplayName', "Modelo reducido")
    
    grid on; box on;
    legend('Interpreter', 'Latex', 'Location', 'Best')
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$v_{rms}$';'[m/s]'},'Interpreter','latex');
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/ApartadoC.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

%%

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
    xlh = get(gca,'xlabel'); xlh.Position(1) = 2000; xlh.Position(2) = 0.0000003;
    
    hold on
    plot(frecuencias, abs(vrms-vrms_CB), ...
        'LineWidth', 1, 'Color', 'k')
  
    grid on; box on;
    xlabel('$Frecuencia$ [Hz]','Interpreter','latex');
    ylabel({'$\Delta v_{rms}$';'[m/s]'},'Interpreter','latex');
    
    if save_flag == 1
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h, './Figures/ApartadoC_error.pdf','-dpdf','-r0','-painters')
    end
    fig = fig+1;
end

