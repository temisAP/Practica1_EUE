classdef barra < material & handle
    %% Propiedades
    properties
        L;
        theta ;
        nodes;
        K;
        T;
        F;
        sigma;
    end
    %% Métodos
     
    methods
        %% Constructor
        function obj = barra(material)
            if nargin > 0
                obj.A = material.A;
                obj.E = material.E;
                obj.alpha = material.alpha;
            end
        end
    end
    
    methods 
        %% Rigidez 
        function val = rigidez(obj)
            obj.K = obj.E * obj.A / obj.L * [ 1  0 -1  0;...
                                              0  0  0  0;...
                                             -1  0  1  0;...
                                              0  0  0  0];
            val = obj.K;
        end
        %% Matriz de rotación
        function T = rotacion(obj)
            alpha = obj.theta;
            T_aux = [cos(alpha), sin(alpha);
                -sin(alpha) cos(alpha)];
            T = [T_aux, zeros(2);
                zeros(2), T_aux];
            obj.T = T;
        end
    end

end