%% Datos 
DT = 50; %K
L = 0.25; 


%% Elementos

% Material
mat = material();
mat.A = j0.0001; %m^2
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
s = estructura(bs);











