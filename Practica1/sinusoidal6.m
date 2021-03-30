%% Datos

% Masa puntual 
M = 50;     %kg

% Viga 
L = 450e-3; %m
b = 20e-3;  %m (perfil cuadrado)

% Excitación
A = 5*9.81; %m/s^2
fmin = 5;   %Hz
fmax = 100; %Hz

% Amortiguamiento
chi = 0.002;

%% Elementos

mat = material();
mat.E = 200e9;  %Pa
mat.rho = 7900; %kg/m^3
mat.nu = 0.3;
