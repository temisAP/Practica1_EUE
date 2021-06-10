function [rms_z, rms_v] = RMS(M, K, F, cond_ini, f)
    % Initial conditions 
    M(cond_ini,:) = []; M(:,cond_ini) = [];
    K(cond_ini,:) = []; K(:,cond_ini) = [];
    F(cond_ini)   = [];

    % Respuesta permanente
    w = 2*pi*f;
    z = F'/(K - w^2*M);
    
    rms_z = sqrt(sum(z.^2)/length(z));
    rms_v = sqrt(sum((w*z).^2)/length(z));

end