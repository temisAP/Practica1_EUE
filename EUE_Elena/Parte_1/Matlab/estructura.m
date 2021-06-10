classdef estructura < handle
    
    %% Propiedades
    
    properties
        gdl;    % Grados de libertad de cada nodo (predeterminado)
        nn;     % Número de nodos
        K;      % Matriz de rigidez
        M;      % Matriz de masa
        F;      % Vector de fuerzas
        N;      % Índices de los nodos que ocupa la estructura
    end
    
    properties (Access = private)
        Matrix; % Para hacer los ensamblajes
        Vector; %
    end
    
    %% Constructor
    
    methods (Access = public)
        
        function obj = estructura(component_s,g)
            
            % Grados de libertad nodales
            obj.gdl = g;
            
            % Descomponer los componentes en elementos nodo a nodo
            elements = descomposicion(component_s);
            
            % Extraer grados de libertad y número de nodos de los
            % componentes
            obj.nn = n_nodes(obj,element_s);            % Número de nodos totales
            GDL = obj.gdl * obj.nn;                     % Grados de libertad totales
            
            
            % Inicializar matrices globales
            if length(size(element_s{1}.M)) == 2
                obj.M = zeros(GDL);
                obj.K = zeros(GDL);
                obj.F = zeros(GDL,1);
            elseif length(size(element_s{1}.M)) == 3
                dimensions = size(element_s{1}.M);
                ww  = dimensions(3);
                obj.M = zeros(GDL,GDL,ww);
            end
            
            % Crear matrices globales
            for e = 1:length(element_s)
                
                % Elemento a ensamblar
                element = element_s{e};
                
                % Ensamblar matriz de masas
                obj.Matrix = element.M;
                obj.M = obj.M + ensamblar(obj,element);
                
                % Ensamblar matriz de rigidez
                obj.Matrix = element.K;
                obj.K = obj.K + ensamblar(obj,element);
                
                % Ensamblar vector de fuerzas
                obj.Vector = element.F;
                obj.F = obj.F + ensamblar_vector(obj,element);
                
            end
        end
    end
    
    %% Otros métodos
    
    methods (Access = private)
        
        % Descomponer los componentes en elementos
        function val = descomposicion(obj,component_s)
            g = obj.gdl;
            for c = 1:length(component_s)
                
                component = component_s{c};
                
                % Esto no vale porque los nodos están cruzados así que hay
                % que hacerlo dejando nodos sin propiedades pero da igual
                % porque luego se van a ensamblar de nuevo quizá
                % interesaría un if para no tener que volver a hacer los
                % cálculos de nuevo ...
                b = reshape(component.M,2*g,2*g,[])
                
                
            end
            
            
        end
        
        % Encontrar el número de nodos totales
        function val = n_nodes(obj,element_s)
            val = 0;
            for b = 1:length(element_s)
                m = max(element_s{b}.nodes);
                if m > val
                    val = m;
                end
            end
        end
        
        
    end
    
    %% Métodos de ensamblaje
    
    methods
        
        % Ensamblaje de matrices
        function [D] = ensamblar(obj,component)
            
            n1 = component.nodes(1);  % nodes: nodos que ocupan a nivel global (vector)
            n2 = component.nodes(2);  % nodes: nodos que ocupan a nivel global (vector)
            g = obj.gdl;              % g    : grados de libertad del componente (escalar)
            nn = obj.nn;              % nn   : número de nodos total (escalar)
            
            
            % Nodos globales que usa cada elemento
            n1 = component.nodes(1);
            n2 = component.nodes(2);
            
            
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
            
            C{n1,n1} = obj.Matrix(A,A);
            C{n2,n1} = obj.Matrix(B,A);
            C{n1,n2} = obj.Matrix(A,B);
            C{n2,n2} = obj.Matrix(B,B);
            
            % Cell 2 matrix
            D = zeros(nn*g);
            D = cell2mat(C);
            
        end
        
        % Ensamblaje de vectores
        function [J] = ensamblar_vector(obj,component)
            
            n1 = component.nodes(1);  % nodes: nodos que ocupan a nivel global (vector)
            n2 = component.nodes(2);  % nodes: nodos que ocupan a nivel global (vector)
            g = obj.gdl;              % g    : grados de libertad locales (escalar)
            nn = obj.nn;              % nn   : número de nodos total (escalar)
            
            
            % Nodos globales que usa cada elemento
            n1 = component.nodes(1);
            n2 = component.nodes(2);
            
            
            % Llenar el cell de ceros
            C = cell(nn,1);
            for i = 1:nn
                C{i} = zeros(g,1);
            end
            
            % Valores que ocupan las submatrices
            A = 1:g;
            B = g+1:2*g;
            
            C{n1} = obj.Vector(A);
            C{n2} = obj.Vector(B);
            
            % Cell 2 matrix
            J = zeros(nn*g,1);
            J = cell2mat(C);
            
        end
    end
    
    %% Métodos de reducción modal
    
    
end