%%% Vigas

classdef viga
    properties
        M           % matriz de masa
        K           % matriz de rigidez
        F           % vector de fuerzas
        n           % nodos
        gdl         % grados de libertad
        frec        % frecuencias naturales
        modos       % modos propios
        psi         % Matriz transformacion CB
    end
end



