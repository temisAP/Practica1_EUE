%% Datos
DT = 50; %K
L = 0.25; %m


%% Elementos

% Material
mat = material();
mat.A = 0.0001; %m^2
mat.E = 70e9; %Pa
mat.alpha = 22.5e-6; %K^-1

% Barras
b12 = barra(mat);
b12.L = L/2^0.5; %m
b12.theta = 0;
b12.nodes = [1 2];

b13 = barra(mat);
b13.L = L; %m
b13.theta = deg2rad(45);
b13.nodes = [1 3];

b23 = barra(mat);
b23.L = L/2^0.5; %m
b23.theta = deg2rad(90);
b23.nodes = [2 3];

%% Estructura

bs = [b12,b13,b23];
s = estructura(bs,2);

%% Desplazamiento térmico

DUt = b13.L * b13.alpha * DT;

a = [1,6];
Ka = s.K(a,a);

v = [cos(b13.theta) sin(b13.theta)]; %Relaciones geométricas

M = [Ka(1,:); v];

Ua = inv(M)*[0;DUt];

%% Reacciones en los apoyos

% Movimientos restringidos
c = [2,3,4,5];

Kca = s.K(c,a);
Fc = Kca * Ua;

%% Tensión sobre cada varilla

U=zeros(length(s.K),1);
U(a) = Ua;
U(c) = zeros(size(c));

% Varilla 1-2

u1 = U(1);
u2 = U(3);

b12.F = b12.K * [u1 0 u2 0]';

% Varilla 1-3

u1 = U(1)*cos(b13.theta) + U(2)*sin(b13.theta);
u2 = U(5)*cos(b13.theta) + U(6)*sin(b13.theta);

b13.F = b13.K * [u1 0 u2 0]';

% Varilla 2-3

u1 = U(4);
u2 = U(6);

b23.F = b23.K * [u1 0 u2 0]';

%% Resultados

disp('Fuerzas de reacción en los apoyos')
disp(Fc)

disp('Desplazamiento de los nodos 1 y 3')
disp(Ua)

disp('Tensión sobre cada varilla')
disp(' b12')
disp(b12.F)
disp(' b13')
disp(b13.F)
disp(' b23')
disp(b23.F)
