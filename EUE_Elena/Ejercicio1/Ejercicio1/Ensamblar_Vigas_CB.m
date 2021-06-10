function [M, K, F, psi] = Ensamblar_Vigas_CB(M1, K1, F1, gdl1, M2, K2, F2, gdl2, psi1, psi2, red1, red2)

    % Inicialization
    n_gdl1 = length(gdl1); % Numero de grados de libertad de viga 1
    n_gdl2 = length(gdl2); % Numero de grados de libertad de viga 2
    M = zeros(n_gdl1 + n_gdl2 - 1);
    K = zeros(n_gdl1 + n_gdl2 - 1);
    psi = zeros(n_gdl1 + n_gdl2 - 1);
    F = zeros(n_gdl1 + n_gdl2 - 1,1);
    
    % Position vector
    pos1 = [1, 3:(n_gdl1+1)]; % Nodo 1 (frontera), nodo 3 interfaz)   
    pos2 = [2:3, (pos1(end)+1):length(M)]; % Nodo 2 (frontera), nodo 3 interfaz) 
    
    % Assemble
    M(pos1,pos1) = M1;
    M(pos2,pos2) = M(pos2,pos2) + M2;
    
    K(pos1,pos1) = K1;
    K(pos2,pos2) = K(pos2,pos2) + K2;
    
    psi(pos1,pos1) = psi1;
    psi(pos2,pos2) = psi(pos2,pos2) + psi2;
    
    F(pos1) = F1;
    F(pos2) = F(pos2) + F2;

    %% REDUCIR
    cond_ini = [1,2];
    M(cond_ini,:) = []; M(:,cond_ini) = [];
    K(cond_ini,:) = []; K(:,cond_ini) = [];
    F(cond_ini) = [];
    psi(cond_ini(1:end),:) = []; psi(:,cond_ini(1:end)) = [];
    
    pos_red1 = n_gdl1-1-red1;
    pos_red2 = length(M)-red2;
    M([pos_red1:(n_gdl1-1) pos_red2:length(M)], :) = []; M(:,[pos_red1:(n_gdl1-1) pos_red2:length(M)]) = []; %Reducir gdl viga 1
    K([pos_red1:(n_gdl1-1) pos_red2:length(K)], :) = []; K(:,[pos_red1:(n_gdl1-1) pos_red2:length(K)]) = [];
    psi(:,[pos_red1:(n_gdl1-1) pos_red2:length(psi)]) = [];
    F([pos_red1:(n_gdl1-1) pos_red2:length(F)]) = [];

    
end