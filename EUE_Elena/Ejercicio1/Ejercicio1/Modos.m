function [modos, frec] = Modos(M,K,cond_ini)
    % Initial conditions 
    M(cond_ini,:) = []; M(:,cond_ini) = [];
    K(cond_ini,:) = []; K(:,cond_ini) = [];
    
    % Calculo frecuencias y modos
    [modos, frec] = eig(K/M);
    frec = sqrt(frec)/(2*pi);
    frec= nonzeros(frec);
    [frec,idx] = sort(frec);
    modos = modos(:,idx);
end