classdef estructura < handle
    %% Propiedades
    properties 
        gdl = 1;    %Grados de libertad de cada nodo (predeterminado)
        K;          %Matriz de rigidez
    end
    %% Métodos 

    
    %% Constructor
    methods (Access = public)
        function obj = estructura(bar_s)
             % Encontrar el número de nodos totales
             nn = 0;
             for b = 1:length(bar_s)
                 m = max(bar_s(b).nodes);
                 if m > nn 
                     nn = m;
                 end
             end
            
            % Construir la matriz de rigidez global
            GDL = obj.gdl * nn;                 % Grados de libertad totales
            obj.K = zeros(GDL);                 % Inicializar matriz de rigidez
            
            for b = 1:length(bar_s)
                bar = bar_s(b);                     % Barra
                K = bar.rigidez();                  % Matriz de rigidez
                T = bar.rotacion();                 % Matriz de rotación
                %bar.K =  T' * bar.K * T;            % Rotación
                obj.K = obj.K + ensamblar(obj,bar,nn);  % Ensamblar matriz
            end
            
        end
    end
    

    
    %% Posicionaje según nodos
    methods
        function [D] = ensamblar(obj,bar,nn)
            K = bar.K;
            n1 = bar.nodes(1);
            n2 = bar.nodes(2);
            g = obj.gdl;
            D=zeros(nn*g);
            if g == 1
                D(n1,n1)=K(1,1);
                D(n1,n2)=K(1,2);
                D(n2,n1)=K(2,1);
                D(n2,n2)=K(2,2);
            elseif g == 2
                disp('2gdls not implemented yet')
            elseif g == 3    
                D(3*n1-2:3*n1,3*n1-2:3*n1)=K(1:3,1:3);
                D(3*n2-2:3*n2,3*n1-2:3*n1)=K(4:6,1:3);
                D(3*n1-2:3*n1,3*n2-2:3*n2)=K(1:3,4:6);
                D(3*n2-2:3*n2,3*n2-2:3*n2)=K(4:6,4:6);
            else
                disp('Not implemented yet')
            end
        end
    end
end