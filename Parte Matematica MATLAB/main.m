clear;
clear all;

%% Condiciones iniciales

gamma = 1.4;
R = 286.9;
theta = 5.352;

u1 = 678;
nu1 = 0;
rho1 = 1.23;
p1 = 1.01e5;
T1 = 286.1;
M1 = 2;


%% Creamos Matrices de Variables

rows = 40;
column = 1;

M = ones(rows,column) * M1;
p = ones(rows,column) * p1;
rho = ones(rows,column) * rho1;
T = ones(rows,column) * T1;
u = ones(rows,column) * u1;

F1 = ones(rows, column) * rho1*u1;
F2 = ones(rows, column) * (rho1*u1*u1 + p1);
F3 = ones(rows, column) * rho1*u1*nu1;
F4 = ones(rows, column) * ((gamma*(gamma-1))*p1*u1 + p1*u1*((u1*u1 + nu1*nu1)/2));

G1 = ones(rows, column) * rho1*u1;
G2 = ones(rows, column) * rho1*u1*nu1;
G3 = ones(rows, column) * (rho1*u1*u1 + p1);
G4 = ones(rows, column) * ((gamma*(gamma-1))*p1*u1 + p1*u1*((u1*u1 + nu1*nu1)/2));

%% Calculos

% Calculamos la primera columna a ver que tal
 