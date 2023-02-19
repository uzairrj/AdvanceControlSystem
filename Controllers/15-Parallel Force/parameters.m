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

q1 = 3;
q2 = -0.2;
q3 = -3;

    H = [[ sin(q1), -cos(q1)*sin(q3), cos(q1)*cos(q3), (2*cos(q1))/5 + (4*cos(q1)*cos(q3))/25 - sin(q1)*(q2 + 3/10)]
[-cos(q1), -sin(q1)*sin(q3), cos(q3)*sin(q1), (2*sin(q1))/5 + cos(q1)*(q2 + 3/10) + (4*cos(q3)*sin(q1))/25]
[       0,         -cos(q3),        -sin(q3),                                              -(4*sin(q3))/25]
[       0,                0,               0,                                                            1]];


eulers = rotm2eul(H(1:3,1:3),"ZYZ");


%operational space desired values
xd = [H(1:3,4); eulers'];

%PD control values
Kp = diag([40,40,40,1,1,1]);
Kd = diag([25,25,25,1,1,1]);
Md = diag([0.12,0.12,0.12,0,0,0]);

%force control parameters
fd = [0,0,3,0,0,0];

Kf = diag([55,55,55,1,1,1]);
Ki = diag([50,50,50,1,1,1]);

%enviroment stiffness
parms.Ke = diag([1,1,1,1,1,1]);

%plane position
parms.plane_pos = [0,0,0.3,0,0,0]';
%parms.plane_pos = [0,0,9.0,0,0,0]';

%plane axis
parms.plane_axis = 3; %z-axis
