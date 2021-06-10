%clear all
close all
fig = 1;
fig_flag = 1;

%% Importar datos de vigas del .pch 

try 
    load('Data.mat')
catch
    ExtractData
end

%% Reducción de Craig-Bampton

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


% Reducción

viga1.CB_reduction();
%viga2.CB_reduction();

%% Ensamblaje

% Renumeración para que encajen
viga2.N = viga2.N + viga1.N(end) - 1;

% Condiciones de contorno

vigas = estructura({viga1,viga2},1);

