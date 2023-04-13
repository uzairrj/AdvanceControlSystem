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

q1 = 1;
q2 = 0.2;
q3 = 0;

    
    H = [[       0,         -cos(q3),        -sin(q3),                                                     -(4*sin(q3))/25]
[ sin(q1), -cos(q1)*sin(q3), cos(q1)*cos(q3),        (2*cos(q1))/5 + (4*cos(q1)*cos(q3))/25 - sin(q1)*(q2 + 3/10)]
[-cos(q1), -sin(q1)*sin(q3), cos(q3)*sin(q1), (2*sin(q1))/5 + cos(q1)*(q2 + 3/10) + (4*cos(q3)*sin(q1))/25 + 3/20]
[       0,                0,               0,                                                                   1]];


eulers = rotm2eul(H(1:3,1:3),"ZYZ");

%operational space desired values
xd = [H(1:3,4); eulers'];

%xd = [0;0;9.00;0;0;0];

%PD control values
Kp = diag([550,550,550,30,30,30]);

Kd = diag([70,70,70,10,10,10]);

%enviroment stiffness
parms.Ke = diag([0,0,850,0,0,0]);

%plane position
parms.plane_pos = [0,0,0.5,0,0,0]';
%parms.plane_pos = [0,0,7.0,0,0,0]';

%plane axis
parms.plane_axis = 3; %z-axis
