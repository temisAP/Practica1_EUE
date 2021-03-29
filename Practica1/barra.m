classdef barra < material
    properties
        L ;
        theta ;
        nodes;
    end
    methods 
        function obj = barra(material)
            if nargin > 0
                obj.A = material.A;
                obj.E = material.E;
                obj.alpha = material.alpha;
            end     
        end
        function K = rigidez
            K = obj.E * obj.A / obj.L * [1 -1;-1 1];
        end
    end 
end