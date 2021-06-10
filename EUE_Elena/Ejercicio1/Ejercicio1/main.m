%%%%%%%%%%% EJERCICIO 1 %%%%%%%%%%%%%%
clear all
clc

%% MODELO COMPLETO
% Viga 1
viga1 = viga();
[viga1.M,viga1.K,viga1.n] = ReadPunchFile_nD('Viga1/viga1.pch',3);
viga1.gdl = 1:length(viga1.n); 
% Viga 2
viga2 = viga();
[viga2.M,viga2.K,viga2.n] = ReadPunchFile_nD('Viga2/viga2.pch',3);
viga2.gdl = 1:length(viga2.n); 
% Ensamblar vigas
vigas = viga();
[vigas.M, vigas.K] = Ensamblar_Vigas(viga1.M,viga1.K, viga1.gdl, ...
                                     viga2.M,viga2.K, viga2.gdl);
vigas.F = zeros(length(vigas.M),1); % Vector de fuerzas externas
vigas.F((viga1.gdl(end)+1)/2)=-1;
% Obtener frecuencias y modos
[viga1.modos, viga1.frec] = Modos(viga1.M, viga1.K, 1); % Modos viga 1
[viga2.modos, viga2.frec] = Modos(viga2.M, viga2.K, 61); % Modos viga 2

cond_ini = [1, length(vigas.M)]; % Nodo 1 y 111 empotrados
[vigas.modos, vigas.frec] = Modos(vigas.M, vigas.K, cond_ini);

% Apartado 1
f = 1:2000; % Rango de frecuencias evaluado [Hz]
for i = 1:length(f)
    [rms_z(i), rms_v(i)] = RMS(vigas.M, vigas.K, vigas.F, cond_ini, f(i));
end

%% CRAIG - BAMPTON
% Viga 1 
viga1_CB = vigas();
nb1 = 1; % Nodo frontera viga 1
nI1 = viga1.gdl(end); % Nodo interfaz
viga1_CB.F= zeros(length(viga1.M),1); % Vector de fuerzas externas
viga1_CB.F((viga1.gdl(end)+1)/2)=-1;
[viga1_CB.M, viga1_CB.K, viga1_CB.F, viga1_CB.psi] = Matrices_CB(viga1.M,...
                                                     viga1.K,viga1_CB.F,nb1, nI1);
% Viga 2 
viga2_CB = vigas();
nI2 = 1; % Nodo interfaz
nb2 = viga2.gdl(end); % Nodo frontera viga 2
viga2_CB.F= zeros(length(viga2.M),1); % Vector de fuerzas externas
[viga2_CB.M, viga2_CB.K, viga2_CB.F, viga2_CB.psi] = Matrices_CB(viga2.M,...
                                                     viga2.K,viga2_CB.F,nb2, nI2);
% Buscar gdl a reducir
[gdl_red1] = Reduccion_CB(viga1.frec);
[gdl_red2] = Reduccion_CB(viga2.frec);

% Ensamblar CB y reducir
vigas_CB = viga();
[vigas_CB.M, vigas_CB.K, vigas_CB.F, vigas_CB.psi] = ...
                                        Ensamblar_Vigas_CB(viga1_CB.M,viga1_CB.K,...
                                        viga1_CB.F,viga1.gdl,viga2_CB.M,viga2_CB.K,...
                                        viga2_CB.F,viga2.gdl,viga1_CB.psi,...
                                        viga2_CB.psi, gdl_red1, gdl_red2);
% Apartado 2
f = 1:2000; % Rango de frecuencias evaluado [Hz]
for i = 1:length(f)
    rms_v_CB(i) = RMS_CB(vigas_CB.M, vigas_CB.K, vigas_CB.F, f(i), vigas_CB.psi);
end

%% FIGURAS
figure(1)
    loglog(f,rms_v)
    hold on
    loglog(f, rms_v_CB)
    grid on
    title('Velocidad RMS')
    legend({'Normal', 'CB'})
    