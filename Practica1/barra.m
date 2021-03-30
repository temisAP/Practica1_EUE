classdef barra < material & handle
    %% Propiedades
    properties
        L;
        theta ;
        nodes;
        K;
        T;
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
        %% Rigidez (método matricial)
        function val = rigidez(obj)
            obj.K = obj.E * obj.A / obj.L * [1 -1;-1 1];
            val = obj.K;
        end
        %% Matriz de rotación
        function T = rotacion(obj)
            alpha = obj.alpha;
            T_aux = [cos(alpha), sin(alpha), 0;
                -sin(alpha) cos(alpha), 0;
                0, 0, 1];
            T = [T_aux, zeros(3);
                zeros(3), T_aux];
            obj.T = T;
        end
    end

end