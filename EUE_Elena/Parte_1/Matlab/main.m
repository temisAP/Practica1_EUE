clear all
close all
fig = 1;
fig_flag = 0;

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
    figure(fig)
        loglog(frecuencias,vrms)
        grid on
        title('vrms')
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

% Estrcuctura completa
vigas_CB = estructura({viga1,viga2},1);

% Calcular las matrices CB a partir de los componentes
vigas_CB.CB_matriz;

% Reducción
f_red = 30;
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
    figure(fig)
        loglog(frecuencias,vrms_CB)
        grid on
        title('vrms')
    fig = fig+1;
end

%% /////// APARTADO C /////// %%


%% Representación gráfica

fig_flag = 1;
if fig_flag == 1
    figure(fig)
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        hold on
        plot(frecuencias, vrms)
        plot(frecuencias, vrms_CB)
        grid on; box on
        title('vrms')
        ylabel('vrms [m/s]')
        xlabel('Frecuencia [Hz]')
        legend('Completa','Reducida')
    fig = fig+1;
end