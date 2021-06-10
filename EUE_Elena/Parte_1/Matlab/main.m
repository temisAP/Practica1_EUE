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

%% Reducción de Craig-Bampton
viga(1).nodes = [1 2];

viga1 = estructura({viga(1)},length(viga(1).N));
%viga1.CB_reduction();

viga(2).nodes = [2 3];
viga2 = estructura({viga(2)},length(viga(2).N));
%viga2.CB_reduction();

%% Ensamblaje

vigas = estructura({viga1,viga2},1);

