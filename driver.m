clear all;
clc;

robot = Robot();

%q1 = 3.3;
%q2 = -0.2;
%q3 = -1.3;

%robot.toBeDeleted();


myrobot = importrobot("RPR_xzx.urdf");

%config = homeConfiguration(myrobot);

%config(1).JointPosition = q1;
%config(2).JointPosition = q2;
%config(3).JointPosition = q3;


%show(myrobot, config);

robot.test(myrobot);

%robot.demo();