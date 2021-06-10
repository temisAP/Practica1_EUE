function [M, K] = Ensamblar_Vigas(M1, K1, gdl1, M2, K2, gdl2)

    % Inicialization
    n_gdl1 = length(gdl1); % Numero de grados de libertad de viga 1
    n_gdl2 = length(gdl2); % Numero de grados de libertad de viga 2
    M = zeros(n_gdl1 + n_gdl2 - 1);
    K = zeros(n_gdl1 + n_gdl2 - 1);
    
    % Position vector
    pos1 = gdl1;
    pos2 = n_gdl1:length(M);
    %pos2(gdl2(1)) = gdl1(end);
    
    % Assemble
    M(pos1,pos1) = M1;
    M(pos2,pos2) = M(pos2,pos2) + M2;
    
    K(pos1,pos1) = K1;
    K(pos2,pos2) = K(pos2,pos2) + K2;

end
