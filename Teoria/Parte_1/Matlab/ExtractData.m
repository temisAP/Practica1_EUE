clc
clear
close all

%% Objetos de las vigas con sus atributos de rigidez, masa y fuerza

% Directorio con los archivos
dir = {'../viga1/viga1.pch','../viga2/viga2.pch'};

% Matrices de masa y rigidez
viga = struct();

for viga_indx = 1:length(dir)
        
    [M,K,N] = ReadPunchFile_nD(dir{viga_indx},[3]);
    
    viga(viga_indx).M = M;
    viga(viga_indx).K = K;
    viga(viga_indx).N = N;
        
end

% Matrices de fuerza
viga(1).F = zeros(length(viga(1).N),1);
viga(1).F(26) = -1;

viga(2).F = zeros(length(viga(2).N),1);

%% Guardar los datos para usarlos luego

save('Data.mat', 'viga')
