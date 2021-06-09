classdef conjunto < handle
    %% Propiedades
    properties 
        gdl;            % Grados de libertad de cada nodo (predeterminado)
        nn;             % Número de nodos
        M;              % Matriz de acoplamientos
        
    end
    %% Constructor
    methods (Access = public)
        function obj = conjunto(component_s,g)
            
            % Extraer grados de libertad y número de nodos de los
            % subconjuntos
            obj.gdl = g;                                % Grados de libertad locales
            obj.nn = n_nodes(obj,component_s);          % Número de nodos totales
            GDL = obj.gdl * obj.nn;                     % Grados de libertad totales
            
            % Inicializar matriz global
            if length(size(component_s{1}.M)) == 2      
                obj.M = zeros(GDL);                     
            elseif length(size(component_s{1}.M)) == 3
                dimensions = size(component_s{1}.M);
                ww  = dimensions(3);
                obj.M = zeros(GDL,GDL,ww);
            end
            
            % Crear matriz global
            for b = 1:length(component_s)
                component = component_s{b};                % Barra
                D = ensamblar(obj,component);
                obj.M = obj.M + ensamblar(obj,component);  % Ensamblar matriz
            end     
        end
    end
    
    %% Otros métodos 
    methods
        % Encontrar el número de nodos totales
        function val = n_nodes(obj,component_s)   
             val = 0;
             for b = 1:length(component_s)
                 m = max(component_s{b}.nodes);
                 if m > val 
                     val = m;
                 end
             end
        end
        % Posicionaje según nodos
        function [D] = ensamblar(obj,component)
            
            n1 = component.nodes(1);  % nodes: nodos que ocupan a nivel global (vector)
            n2 = component.nodes(2);  % nodes: nodos que ocupan a nivel global (vector)
            g = obj.gdl;              % g    : grados de libertad locales (escalar)
            nn = obj.nn;              % nn   : número de nodos total (escalar)
            

            % Nodos globales que usa cada elemento
            n1 = component.nodes(1);
            n2 = component.nodes(2);
            
            if length(size(component.M)) == 2         
                % Llenar el cell de ceros
                C = cell(nn);
                for i = 1:nn
                    for j = 1:nn
                       C{i,j} = zeros(g,g);
                     end
                end

                % Valores que ocupan las submatrices
                A = 1:g;
                B = g+1:2*g;

                C{n1,n1} = component.M(A,A);
                C{n2,n1} = component.M(B,A);
                C{n1,n2} = component.M(A,B);
                C{n2,n2} = component.M(B,B);  

                % Cell 2 matrix
                D = zeros(nn*g);
                D = cell2mat(C);
                
            elseif length(size(component.M)) == 3   % Un tensor que depende de la frecuencia o el tiempo
                
                dimensions = size(component.M);
                ww  = dimensions(3);
                
                % Llenar el cell de ceros
                C = cell(nn,nn,ww);
                for i = 1:nn
                    for j = 1:nn
                        for k = 1:ww
                            C{i,j,k} = zeros(g,g);
                        end
                     end
                end

                % Valores que ocupan las submatrices
                A = 1:g;
                B = g+1:2*g;
                
                for k = 1:ww
                    C{n1,n1,k} = component.M(A,A,k);
                    C{n2,n1,k} = component.M(B,A,k);
                    C{n1,n2,k} = component.M(A,B,k);
                    C{n2,n2,k} = component.M(B,B,k);  
                end

                % Cell 2 matrix
                D = zeros(nn*g,nn*g,ww);
                D = cell2mat(C);
            end
                
        end
        function [J] = ensamblar_vector(obj,component)
            
            n1 = component.nodes(1);  % nodes: nodos que ocupan a nivel global (vector)
            n2 = component.nodes(2);  % nodes: nodos que ocupan a nivel global (vector)
            g = obj.gdl;              % g    : grados de libertad locales (escalar)
            nn = obj.nn;              % nn   : número de nodos total (escalar)
            

            % Nodos globales que usa cada elemento
            n1 = component.nodes(1);
            n2 = component.nodes(2);
            
            if length(size(component.V)) == 2         
                % Llenar el cell de ceros
                C = cell(nn,1);
                for i = 1:nn
                    C{i} = zeros(g,1);
                end

                % Valores que ocupan las submatrices
                A = 1:g;
                B = g+1:2*g;

                C{n1} = component.V(A);
                C{n2} = component.V(B);

                % Cell 2 matrix
                J = zeros(nn*g,1);
                J = cell2mat(C);
                
            elseif length(size(component.V)) == 3   % Un tensor que depende de la frecuencia o el tiempo
                
                dimensions = size(component.V);
                ww  = dimensions(3);
                
                % Llenar el cell de ceros
                C = cell(nn,ww);
                for i = 1:nn
                    for j = 1:ww
                        C{i,j} = zeros(g,ww);
                    end
                end

                % Valores que ocupan las submatrices
                A = 1:g;
                B = g+1:2*g;
                
                for k = 1:ww
                    C{n1,k} = component.V(A,k);
                    C{n2,k} = component.V(B,k); 
                end

                % Cell 2 matrix
                J = zeros(nn*g,ww);
                J = cell2mat(C);
            end
                
        end
    end
end