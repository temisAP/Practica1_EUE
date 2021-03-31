classdef estructura < handle
    %% Propiedades
    properties 
        gdl;    % Grados de libertad de cada nodo (predeterminado)
        nn;     % Número de nodos
        K;      % Matriz de rigidez
    end
    %% Constructor
    methods (Access = public)
        function obj = estructura(bar_s,g)
            obj.gdl = g;
            obj.nn = obj.n_nodes(bar_s);
            GDL = obj.gdl * obj.nn;                 % Grados de libertad totales
            obj.K = zeros(GDL);                     % Inicializar matriz de rigidez
            
            for b = 1:length(bar_s)
                bar = bar_s(b);                      % Barra
                bar.rigidez();                       % Matriz de rigidez
                bar.rotacion();                      % Matriz de rotación
                obj.K = obj.K + ensamblar(obj,bar);  % Ensamblar matriz
            end     
        end
    end
    
    %% Otros métodos 
    methods
        % Encontrar el número de nodos totales
        function val = n_nodes(obj,bar_s)   
             val = 0;
             for b = 1:length(bar_s)
                 m = max(bar_s(b).nodes);
                 if m > val 
                     val = m;
                 end
             end
        end
        % Posicionaje según nodos
        function [D] = ensamblar(obj,bar)
            K = bar.T' * bar.K * bar.T; %Rotada
            n1 = bar.nodes(1);
            n2 = bar.nodes(2);
            g = obj.gdl;
            nn = obj.nn;
            
            C = cell(nn);
            for i = 1:nn
                for j = 1:nn
                   C{i,j} = zeros(g,g);
                 end
            end
            
            D=zeros(nn*g);
            if g == 1
                C{n1,n1} = K(1,1);
                C{n2,n1} = K(2,1);
                C{n1,n2} = K(1,2);
                C{n2,n2} = K(2,2); 
            elseif g == 2
                C{n1,n1} = K(1:2,1:2);
                C{n2,n1} = K(3:4,1:2);
                C{n1,n2} = K(1:2,3:4);
                C{n2,n2} = K(3:4,3:4);    
            elseif g == 3    
                C{n1,n1} = K(1:3,1:3);
                C{n2,n1} = K(4:6,1:3);
                C{n1,n2} = K(1:3,4:6);
                C{n2,n2} = K(4:6,4:6);  
            else
                disp('Not implemented yet')
            end
            D = cell2mat(C);
        end
    end
end