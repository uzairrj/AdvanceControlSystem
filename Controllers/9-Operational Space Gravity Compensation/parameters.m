clear;
clc;

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


%q1 = 2;
%q2 = -0.1;
%q3 = -3;

%singularity
q1 = -2;
q2 = 0.12;
q3 = 3;

    H = [[ sin(q1), -cos(q1)*sin(q3), cos(q1)*cos(q3), (2*cos(q1))/5 + (4*cos(q1)*cos(q3))/25 - sin(q1)*(q2 + 3/10)]
[-cos(q1), -sin(q1)*sin(q3), cos(q3)*sin(q1), (2*sin(q1))/5 + cos(q1)*(q2 + 3/10) + (4*cos(q3)*sin(q1))/25]
[       0,         -cos(q3),        -sin(q3),                                              -(4*sin(q3))/25]
[       0,                0,               0,                                                            1]];


eulers = rotm2eul(H(1:3,1:3),"ZYZ");


%operational space desired values
xd = [H(1:3,4); eulers'];

%PD control values
Kp = diag([1000,1000,1000,5,5,5]);

Kd = diag([250,250,250,2,2,2]);

