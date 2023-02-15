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

%Joint variables values
qd = [2.1, -0.2, 1];

%PD control values
Kp = [[150,0,0];
      [0,150,0];
      [0,0,150]];

Kd = [[50,0,0];
      [0,50,0];
      [0,0,50]];