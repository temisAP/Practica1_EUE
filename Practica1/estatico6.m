%% Datos 
DT = 50; %K
L = 0.25; 


%% Objetos

% Material
mat = material();
mat.A = 0.0001; %m^2
mat.E = 70e9; %Pa
mat.alpha = 22.5e-6; %K^-1

% Barras
b12 = barra(mat);
b12.L = 0.25/2^0.5;
b12.theta = 0;
b12.nodes = [1 2];

b13 = barra(mat);
b13.L = 0.25;
b13.theta = deg2rad(45);
b13.nodes = [1 3];

b23 = barra(mat);
b23.L = 0.25/2^0.5;
b23.theta = deg2rad(90);
b23.nodes = [2 3];

%% Ensamblaje

bs = [b12,b13,b23];











