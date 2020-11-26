clear;
clear all;

%% Condiciones iniciales

gamma = 1.4;
R = 286.9;
theta = 5.352;

u1 = 678;
v1 = 0;
rho1 = 1.23;
p1 = 1.01e5;
T1 = 286.1;
M1 = 2;

C = 0.5;
E = 10;
x = 0;
H = 40;



%% Creamos Matrices de Variables en la primera columna

rows = 40;
column = 1;

M = ones(rows,column) * M1;
p = ones(rows,column) * p1;
rho = ones(rows,column) * rho1;
T = ones(rows,column) * T1;
u = ones(rows,column) * u1;

F1 = ones(rows, column) * rho1*u1;
F2 = ones(rows, column) * (rho1*u1*u1 + p1);
F3 = ones(rows, column) * rho1*u1*v1;
F4 = ones(rows, column) * ((gamma*(gamma-1))*p1*u1 + p1*u1*((u1*u1 + v1*v1)/2));

G1 = ones(rows, column) * rho1*u1;
G2 = ones(rows, column) * rho1*u1*v1;
G3 = ones(rows, column) * (rho1*u1*u1 + p1);
G4 = ones(rows, column) * ((gamma*(gamma-1))*p1*u1 + p1*u1*((u1*u1 + v1*v1)/2));


%% Calculos

deltaY = H / rows;


j = 1;

while(x<=65)
    
    i = 1;
    for i = 1:40
        
    % Para calcular p, rho, u,v,etc... per la columna 2 perimero nos vamos
    % a la columna 1
    
    % 1 Con theta, M, etc... calculamos mu,
    mu = arcsind(1 / M(j,i));
    deltaX = deltaY / Max(Abs(theta - mu), Abs(thata + mu));
    deltaETA = C * deltaX;
    
    % 2 Con mu, deltaX, deltaETA, etc... calculamos X y con x calculamos h y
    % ys
    x = x + deltaETA;
    h = h(x, E, theta, H);
    y_s = y_s(x,E,theta);
    
    % 3 Con h y y_s calculamos ETA y NU
    ETA = x;
    NU = (deltaY*j - y_s) / h;
    
    % 4 con NU, h calculamos parcial NU / parcial X
    parcialNU_parcialX = parcial_NU_parcial_x(NU, theta, h, E);
    
    % 5 Con parcialNU_parcialX calculamos parcialF_parcialETA para F1234
    
    
        
        
    end

end