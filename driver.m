clear all

robot = Robot();

myrobot = importrobot("RPR_xzx.urdf");

show(myrobot)

y = robot.test(0.02, 0.1256, 0.4, 0.2)