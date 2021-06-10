function [M_CB, K_CB, F_CB, psi] = Matrices_CB(M,K,F,nb,nI)
    % Vector de nodos interiores
    n_int = 2:(length(M)-1);
    n = [nb, nI, n_int];
    % Ordenar matrices
    M = M(n,n);
    Mff = M(1:2,1:2); % Masa nodos frontera
    Mii = M(3:end,3:end); % Masa nodos interior
    K = K(n,n);
    Kif = K(3:end,1:2);
    Kii = K(3:end,3:end);
    F = F(n);
    % Transformar a espacio Craig-Bampton
        % Calculo modos estaticos
        phi_s = [eye(length(Mff));...
                 -inv(Kii)*Kif];
        % Calculo modos elasticos interiores
        [mod, frec] = Modos(Mii, Kii,[]); % Al ser nodos interiores no hay ci
        phi_i = [zeros(length(Mff),length(mod)); ...
                 mod];
        % Calculo matriz de transformacion     
        psi = [phi_s phi_i];
        % Matrices Craig-Bampton
        M_CB = psi'*M*psi;
        K_CB = psi'*K*psi;
        F_CB = psi'*F;
end
