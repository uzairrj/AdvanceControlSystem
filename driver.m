clear all;

robot = Robot();

myrobot = importrobot("RPR_xzx.urdf");

%show(myrobot);

robot.test(myrobot);

%robot.demo();