%sine parameters
A = 1;
Fc = 1;

%lengths
lenghts_values = [[0.02, 0.1256, 0.4]; 
                  [0.03, 0.03, 0.3]; 
                  [0.02, 0.1256, 0.16]];
density = 2700; %Almunium density
mass_values = [pi * lenghts_values(1,1)^2*lenghts_values(1,3)*density, lenghts_values(2,1)*lenghts_values(2,2)*lenghts_values(3,3)*density, lenghts_values(3,1)^2*lenghts_values(3,3)*density];
           

%masses
params.m = mass_values;

%gravity
params.g = 9.806;

%Friction
params.Fv = [0,0,0];

%PD control values
%Kp = diag([90,90,90,0,0,0]);

%Kd = diag([10,10,10,0,0,0]);

%Md = diag([0.2,0.2,0.2,0,0,0]);


Kp = diag([790,790,790,70,70,70]);

Kd = diag([170,170,170,50,50,50]);

Md = diag([0.12,0.12,0.12,0.01,0.01,0.01]);


Kpt = diag([350,350,350,50,50,50]);
Kdt = diag([125,125,125,15,15,15]);
Mt = diag([50,50,50,10,10,10]);


%enviroment stiffness
parms.Ke = diag([0,0,190,0,0,0]);

%plane position
parms.plane_pos = [0,0,0.4,0,0,0]';
%parms.plane_pos = [0,0,9.0,0,0,0]';

%plane axis
parms.plane_axis = 3; %z-axis
