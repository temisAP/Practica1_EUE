function [rms_v] = RMS_CB(M, K, F, f, psi)

    % Respuesta permanente
    w = 2*pi*f;
    z = psi*((K - w^2*M)\F);
    
    rms_z = sqrt(sum(z.^2)/length(z));
    rms_v = sqrt(sum((w*z).^2)/length(z));

end