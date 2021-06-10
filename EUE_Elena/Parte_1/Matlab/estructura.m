classdef estructura < handle
    
    %% Propiedades
    
    properties
        
        % Nodos
        
        gdl;    % Grados de libertad de cada nodo (predeterminado)
        GDL;    % Grados de libertad totales
        NN;     % Número de nodos totales de la estructura
        N;      % Índices de los nodos que ocupa la estructura
        
        bn      % Nodos de frontera (boundary nodes)
        in      % Nodos internos    (inner nodes)
        iface   % Nodos de interfaz con otros elementos
        
        % Propiedades de la estructura
        
        component_s; % Componentes de la estructura
        
        K;      % Matriz de rigidez
        M;      % Matriz de masa
        F;      % Vector de fuerzas
                
        eigval; % Autovalores
        eigvec; % Autovectores
        
        q;      % Vector de estado
        
        % Para la reducción CB
        
        K_CB;   % Matriz de rigidez de Craig-Bampton
        M_CB;   % Matriz de masa de Craig-Bampton
        F_CB;   % Vector de fuerzas de Craig-Bampton
        psi;    % Matriz modal
        gdl_red;% Grados de libertad a reducir según criterio
        
    end
    
    properties (Access = private)
        
        % Para hacer los ensamblajes
        Matrix; 
        Vector;
        
        % Bandas de frecuencia
        Lf = [14.1;17.8;22.4;28.2;35.5;44.7;56.2;70.8;89.1;112;141;178;224;...
            282;355;447;562;708;891;1122;1413;1778;2239;2818;3548;4467;5623;7079;...
            8913;11220;14130;17780];
        Uf = [17.8;22.4;28.2;35.5;44.7;56.2;70.8;89.1;112;141;178;224;282;355;...
            447;562;708;891;1122;1413;1778;2239;2818;3548;4467;5623;7079;...
            8913;11220;14130;17780;22390];
        
    end
    
    %% Constructor
    
    methods (Access = public)
        
        function obj = estructura(component_s,g)
            
            disp('Creando estructura...')
            
            % Compoenentes de la estructura
            obj.component_s = component_s;
            
            % Grados de libertad nodales
            obj.gdl = g;
            
            % Extraer grados de libertad y número de nodos de los componentes
            obj.NN  = N_nodes(obj);             % Número de nodos totales
            obj.GDL = obj.gdl * obj.NN;                     % Grados de libertad totales
            disp(['    Número de nodos totales: ',num2str(obj.NN)]);
            disp(['    Número de grados de libertad totales: ',num2str(obj.GDL)]);
            
            % Nodos de la estructura
            disp('    Asignando nodos');
            asignar_nodos(obj);
            
            % Crear matrices globales
            disp('    Ensamblando matrices');
            matrices_globales(obj);
                                    
            disp('Estructura creada')
            disp('')
            
        end
    end
    
    %% Otros métodos
    
    methods (Access = private)
        
        % Encontrar el número de nodos totales de la estrcutura
        function val = N_nodes(obj)
            
            component_s = obj.component_s;           
            
            val = 0;
            for b = 1:length(component_s)
                mx = max(component_s{b}.N);
                if mx > val
                    val = mx;
                end
            end
        end
        
        % Asignar los nodos a la estructura
        function asignar_nodos(obj)
            
            component_s = obj.component_s;
            
            for c = 1:length(component_s)
                % Nodos totales
                obj.N  = [obj.N ; component_s{c}.N];
                % Nodos frontera
                try
                    obj.bn = [obj.bn ; component_s{c}.bn];
                catch
                    disp('    *No se han encontrado nodos frontera')
                end
                % Nodos interiores
                try
                    obj.in = [obj.in ; component_s{c}.in];
                catch
                    disp('    *No se han encontrado nodos interiores')
                end
                % Nodos interfaz
                try
                    obj.iface = [obj.iface ; component_s{c}.iface];
                catch
                    disp('    *No se han encontrado nodos de interfaz')
                end
            end
            
            obj.N = unique(obj.N);
            
            % Autocompletar
            if isempty(obj.in) && ~isempty(obj.bn)
                obj.in = obj.N(obj.N~=obj.bn & obj.N~=obj.iface);
            elseif ~isempty(obj.in) && ~isempty(obj.bn)
                for c = 1:length(component_s)
                    obj.bn = [obj.bn; component_s{c}.bn];
                    obj.in = [obj.in; component_s{c}.in];
                end
            end
            
            obj.in = unique(obj.in);
            obj.bn = unique(obj.bn);
            obj.iface = unique(obj.iface);
            
        end
        
        % Crear matrices globales a partir de los componentes
        function matrices_globales(obj)
            
            component_s = obj.component_s;
            
            GDL = obj.GDL;
            
            % Inicializar matrices globales
            obj.M = zeros(GDL);
            obj.K = zeros(GDL);
            obj.F = zeros(GDL,1);
            
            for c = 1:length(component_s)
                
                % Elemento a ensamblar
                component = component_s{c};
                
                % Ensamblar matriz de masas
                obj.Matrix = component.M;    
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
        
        function [eigval, eigvec] = Modos(obj,M,K)
            
            
            % Si no se especifica qué matrices son son las de la estructura
            % con sus condiciones de contorno
            if nargin == 1
                disp('    Modos de la estructura')
                M = obj.M;
                K = obj.K;
                
                cc = [obj.bn; obj.iface];
                
                % Condiciones de contorno
                M(cc,:) = []; M(:,cc) = [];
                K(cc,:) = []; K(:,cc) = [];
                
            end
            
            % Autovalores y autovectores
            [eigvec, eigenval_matrix] = eig( M\K );
            eigenval_matrix = sqrt(eigenval_matrix);
            
            eigenval_matrix = eigenval_matrix(sub2ind(size(eigenval_matrix),...
                1:size(eigenval_matrix,1),1:size(eigenval_matrix,2)));
            
            [eigval,idx] = sort(eigenval_matrix);
            
            eigval = eigval/(2*pi);
            eigvec = eigvec(:,idx);
            
            if nargin == 1
                obj.eigval = eigval;
                obj.eigvec = eigvec;
            end
            
        end
            
    function q = dinamic_solution(obj,f,M,K,F)
            
            w = 2*pi*f; % rad/s
            
            % Si no se especifica qué matrices son son las de la estructura
            % con sus condiciones de contorno
            if nargin == 2
                M = obj.M;
                K = obj.K;
                F = obj.F;
             
                % Condiciones de contorno
                cc = sort([obj.bn; obj.iface]);
                aa = setdiff(obj.N,cc);
                
                M(cc,:) = []; M(:,cc) = [];
                K(cc,:) = []; K(:,cc) = [];
                F(cc)   = [];
                               
            end
            
            % Respuesta permanente
            z = inv(K - w^2*M)*F;
            q = zeros(obj.GDL,1);
            q(aa) = z;
            obj.q = q;
        
    end
    
    end
    
    
    %% Métodos de reducción modal
    methods
        
        function CB(obj)
            CB_transformation(obj);
            CB_getGDLred(obj);
            CB_reduction(obj);          
        end
        
        function CB_transformation(obj)
            
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

            M = CB_assembly(obj,M);
            K = CB_assembly(obj,K);
                            
            Mff = cell2mat({M{1,1},M{1,2};M{2,1},M{2,2}});
            Mii = M{3,3};                       
            
            Kif = cell2mat({K{3,1},K{3,2}});
            Kii = K{3,3};
            
            F = F(gdl_CB);
            
            % Matriz de transformación CB
            phi_s = [eye(size(Mff));...
                -inv(Kii)*Kif];
            
            [~, modos] = Modos(obj,Mii,Kii);
            
            phi_i = [zeros(length(Mff),size(modos,1)); modos];
            
            psi = [phi_s, phi_i];
            
            % Convertir a matriz
            M = cell2mat(M);
            K = cell2mat(K);
            
            % Return
            obj.M_CB = psi'*M*psi;
            obj.K_CB = psi'*K*psi;
            obj.F_CB = psi'*F;
            obj.psi  = psi;
            
        end
               
        
        function CB_getGDLred(obj,f_red)

            % Grados de libertad a reducir de cada componente
            
            Lf = obj.Lf;                % Límites inferiores de bandas de frecuencia
            Uf = obj.Uf;                % Límites superiores de bandas de frecuencia
            k = zeros(length(Lf),1);
            
            try
                eigval = obj.eigval;
            catch
                Modos(obj);
                eigval = obj.eigval;
            end

            % Modos por banda
            for i = 1: length(k)
                k(i) = length(eigval(eigval>Banda_baja(i) & eigval<Banda_alta(i)));
            end
            % posiciones con más de f_red modos
            pos = find(k>f_red);
            % Nos quedamos con solo esas posiciones en las matrices
            gdl_red = length((find(eigval>=Banda_baja(pos(1)))));
            obj.gld_red = gdl_red;
                                 
        end
        
        % Crear matrices reducidas globales a partir de los componentes
        %       Si no se llama esta función se usan directamente las 
        %       matrices reducidas de la estructura, si las hubiera.
        function CB_matriz(obj)
            
            component_s = obj.component_s;
            
            for c = 1:length(component_s)
                gdl(c) = component_s{c}.GDL;
            end
             
            try
                
                % Matrices de las dos vigas
                K1 = component_s{1}.K_CB;
                M1 = component_s{1}.M_CB;
                F1 = component_s{1}.F_CB;
                psi_1 = component_s{1}.psi;
                
                K2 = component_s{2}.K_CB;
                M2 = component_s{2}.M_CB;
                F2 = component_s{2}.F_CB;
                psi_2 = component_s{2}.psi;
                

                % Matrices globales
                M = zeros(gdl(1) + gdl(2) - 1);
                K = zeros(gdl(1) + gdl(2) - 1);
                F = zeros(gdl(1) + gdl(2) - 1, 1);
                psi = zeros(gdl(1) + gdl(2) - 1);

                % Vector de posiciones
                pos_1 = [1, 3:(gdl(1)+1)];
                pos_2 = [2:3, (pos_1(end)+1):length(M)];

                % Ensamblar matrices
                M(pos_1,pos_1) = M(pos_1,pos_1) + M1;
                M(pos_2,pos_2) = M(pos_2,pos_2) + M2;

                K(pos_1,pos_1) = K(pos_1,pos_1) + K1;
                K(pos_2,pos_2) = K(pos_2,pos_2) + K2;

                F(pos_1) = F(pos_1) + F1;
                F(pos_2) = F(pos_2) + F2;

                psi(pos_1,pos_1) = psi(pos_1,pos_1) + psi_1;
                psi(pos_2,pos_2) = psi(pos_2,pos_2) + psi_2;

                obj.M_CB = M;
                obj.K_CB = K;
                obj.F_CB = F;
                obj.psi  = psi;
            catch
                disp('De momento esta parte solo soporta dos vigas con un nodo de interfaz')
            end
            
        end
        
        function CB_reduction(obj)
            
            component_s = obj.component_s;
            
            nodos_frontera = [1 2];
            nodos_interfaz = [3];
            
            for c = 1:length(component_s)
                gdl_red(c) = component_s{c}.gdl_red;
            end
            
        end        
    end
    
    methods (Access = private)
        function [C] = CB_assembly(obj,M,gdl_CB)
            
            if nargin == 2
                in = obj.in;
                bn = obj.bn;
                iface = obj.iface;
            else
                if isempty(gdl_B)
                    in = length(obj.in);
                    bn = length(obj.bn);
                    iface = length(obj.iface);         
                else
                    in = gdl_CB{1};
                    bn = gdl_CB(2);
                    iface = gdl_CB(3);
                end
            end
                        
            C = cell(3);
            C{1,1} = M(bn,bn);
            C{1,2} = M(bn,iface);
            C{1,3} = M(bn,in);
            
            C{2,1} = M(iface,bn);
            C{2,2} = M(iface,iface);
            C{2,3} = M(iface,in);
            
            C{3,1} = M(in,bn);
            C{3,2} = M(in,iface);
            C{3,3} = M(in,in);
                        
        end

    end
    
end