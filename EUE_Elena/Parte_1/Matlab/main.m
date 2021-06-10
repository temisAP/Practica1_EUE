clear all
close all
fig = 1;
fig_flag = 1;

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

vigas = estructura({viga1,viga2},1);

%% Velocidad rms 

frecuencias = linspace(1,2000,2000); % Hz

for f=1:length(frecuencias)
    q = vigas.dinamic_solution(frecuencias(f));
    qrms(f) = rms(q);
end

vrms = qrms * 2 * pi .*frecuencias;

%% /////// APARTADO B /////// %%
%% Reducción de la estructura

% Transformación Craig-Bampton de las vigas
viga1.CB_transformation();
viga2.CB_transformation();

% Estrcuctura completa
vigas_CB = estructura({viga1,viga2},1);

% Reducción
f_red = 30;
vigas_CB.reduccion(f_red,{viga1,viga2});

%% Velocidad rms

frecuencias = linspace(1,2000,2000); % Hz

for f=1:length(frecuencias)
    q = vigasCB.dinamic_solution(frecuencias(f));
    qrms(f) = rms(q);
end

vrms = qrms * 2 * pi .*frecuencias;
