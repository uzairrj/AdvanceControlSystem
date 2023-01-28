%sine parameters
A = 1;
Fc = 1;

%Joint variables values
qd = [1.3, -0.5, -1.3];

%PD control values
Kp = [[100,0,0];
      [0,100,0];
      [0,0,100]];

Kd = [[10,0,0];
      [0,10,0];
      [0,0,10]];