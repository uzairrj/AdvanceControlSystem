%mass
density = 2.7;
La = 1.1;
Lb = 1.1;
Lc = 3.5;
mass = La*Lb*Lc*density;

%mass = 1;

%gravity
g = 9.806;

%mid point
d = Lc/2;

%d = 0.2;

%1 DoF dynamics
params.I = 0.1;
params.F = 0.2;
params.G = mass * g * d;

%kpi
K = 30.2;

%lambda
lambda = 100;

%PD control values
Kp = 25; %it is not used

Kd = 500; 