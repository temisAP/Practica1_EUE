classdef estructura < handle
    
    %% Propiedades
    
    properties
        gdl;    % Grados de libertad de cada nodo (predeterminado)
        NN;     % Número de nodos totales de la estructura
        N;      % Índices de los nodos que ocupa la estructura
        
        bn      % Nodos de frontera (boundary nodes)
        in      % Nodos internos    (inner nodes)
        iface   % Nodos de interfaz con otros elementos
        
        K;      % Matriz de rigidez
        M;      % Matriz de masa
        F;      % Vector de fuerzas
        
        eigval; % Autovalores
        eigvec; % Autovectores
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
            
            % Extraer grados de libertad y número de nodos de los componentes
            obj.NN = N_nodes(obj,component_s);          % Número de nodos totales
            disp(['Número de nodos de la estructura: ',num2str(obj.NN)]);
            GDL = obj.gdl * obj.NN;                     % Grados de libertad totales
            
            % Nodos de la estructura
            for c = 1:length(component_s)
                % Nodos totales
                obj.N  = [obj.N ; component_s{c}.N];
                % Nodos frontera
                try 
                    obj.bn = [obj.bn ; component_s{c}.bn];
                catch
                    disp('No se han encontrado nodos frontera')
                end
                % Nodos interiores
                try
                    obj.in = [obj.in ; component_s{c}.in];
                catch
                    disp('No se han encontrado nodos interiores')
                end
                % Nodos interfaz
                try
                    obj.iface = [obj.iface ; component_s{c}.iface];
                catch
                    disp('No se han encontrado nodos de interfaz')
                end
            end
            
            % Autocompletar 
            if isempty(obj.in) && ~isempty(obj.bn)
                obj.in = obj.N(obj.N~=obj.bn & obj.N~=obj.iface);
            elseif ~isempty(obj.in) && ~isempty(obj.bn)
                obj.bn = [obj.bn; component_s{c}.bn];
                obj.in = [obj.in; component_s{c}.in];
            end
            
            % Inicializar matrices globales
            obj.M = zeros(GDL);
            obj.K = zeros(GDL);
            obj.F = zeros(GDL,1);
            
            % Crear matrices globales
            for c = 1:length(component_s)
                
                % Elemento a ensamblar
                component = component_s{c};
                
                % Ensamblar matriz de masas
                obj.Matrix = component.M;
                D = ensamblar(obj,component);
                obj.M = obj.M + ensamblar(obj,component);
                
                % Ensamblar matriz de rigidez
                obj.Matrix = component.K;
                obj.K = obj.K + ensamblar(obj,component);
                
                % Ensamblar vector de fuerzas
                obj.Vector = component.F;
                obj.F = obj.F + ensamblar_vector(obj,component);
                
            end
        end
    end
    
    %% Otros métodos
    
    methods (Access = private)
        
        % Encontrar el número de nodos totales de la estrcutura
        function val = N_nodes(obj,component_s)
            val = 0;
            for b = 1:length(component_s)
                mx = max(component_s{b}.N);
                if mx > val
                    val = mx;
                end
            end
        end
        
        
    end
    
    %% Métodos de ensamblaje
    
    methods (Access = private)
        
        % Ensamblaje de matrices
        function [D] = ensamblar(obj,component)
            
            g = obj.gdl;                % g    : grados de libertad del componente (escalar)
            NN = obj.NN;                % nn   : número de nodos total (escalar)
            
            % Nodos globales que usa cada elemento
            n = component.N(:);
            
            % Separar obj.Matrix en submatrices por nodo
            sub_matrix_number = length(obj.Matrix)/g;           % Número de submatrices
            sub_matrix_sizes  = ones(sub_matrix_number,1) * g;  % Tamaño de las submatrices
            
            M = mat2cell(obj.Matrix,sub_matrix_sizes,sub_matrix_sizes);
            
            % Llenar el cell de ceros
            C = cell(NN);
            for i = 1:NN
                for j = 1:NN
                    C{i,j} = zeros(g,g);
                end
            end
            
            % Valores que ocupan las submatrices
            
            for i = 1:sub_matrix_number
                for j = 1:sub_matrix_number
                    C{n(i),n(j)} = M{i,j};
                end
            end
            
            % Cell 2 matrix
            D = cell2mat(C);
            
        end
        
        % Ensamblaje de vectores
        function [J] = ensamblar_vector(obj,component)
            
            g = obj.gdl;              % g    : grados de libertad locales (escalar)
            NN = obj.NN;              % nn   : número de nodos total (escalar)
            
            
            % Nodos globales que usa cada elemento
            n = component.N(:);
            
            % Separar obj.Matrix en submatrices por nodo
            
            sub_vector_number = length(obj.Matrix)/g;           % Número de submatrices
            sub_vector_size1  = ones(sub_vector_number,1) * g;  % Tamaño de las submatrices
            
            V = mat2cell(obj.Vector,sub_vector_size1);
            
            
            % Llenar el cell de ceros
            C = cell(NN);
            for i = 1:NN
                C{i} = zeros(g,1);
            end
            
            % Valores que ocupan las submatrices
            for i = 1:sub_vector_number
                C{n(i)} = V{i};
            end
            
            
            % Cell 2 matrix
            J = cell2mat(C);
            
        end
    end
    
    %% Métodos de análisis dinámico
    
    methods
        
        function [eigenval, eigenvect] = modos(obj)
            
            eigenval = 0;
            eigenvect = 0;
            
        end
        
    end
    
    
    %% Métodos de reducción modal
    methods
        
        function CB_reduction(obj)
            
            % Nodos interiores
            in = obj.in';
            nb = obj.bn';
            iface = obj.iface';
            
            % Matrices de la estructura
            M = obj.M;
            K = obj.K;
            F = obj.F;
            
            % Grados de libertad
            gdl_CB = [nb, iface, in];
            
            % Matrices de masa y rigidez ordenadas
            M = M(gdl_CB,gdl_CB);
            Mff = M(1:2,1:2);
            Mii = M(3:end,3:end);
            
            K = K(gdl_CB,gdl_CB);
            Kif = K(3:end,1:2);
            Kii = K(3:end,3:end);
            
            F = F(gdl_CB);
            
            % Matriz de transformación CB
            phi_s = [eye(size(Mff));...
                -inv(Kii)*Kif];
            
            % Obtencion de los modos y frecuencias propias^2
            [mod_prop, frec_matrix] = eig( M\K );
            frec_matrix = frec_matrix(sub2ind(size(frec_matrix),...
                1:size(frec_matrix,1),1:size(frec_matrix,2)));
            
            [frec_prop,idx] = sort(frec_matrix);
            frec_prop = frec_prop/(2*pi);
            mod_prop = mod_prop(:,idx);
            
            phi_i = [zeros(length(Mff),size(mod_prop,1)); mod_prop];
            
            psi = [phi_s, phi_i];
            
            % Return
            obj.M = psi'*M*psi;
            obj.K = psi'*K*psi;
            obj.F = psi'*F;
            
        end
        
    end
end