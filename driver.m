clear all

robot = Robot();

myrobot = importrobot("RPR_xzx.urdf");

config = homeConfiguration(myrobot);


t1 = 3*(pi/180);
t2 = 1*(pi/180);
config(1).JointPosition = t1
config(2).JointPosition = 6
config(3).JointPosition = t2

jacob = geometricJacobian(myrobot,config,"ee")

%robot.TransoformationMatrix

robot.solveAnalyticalJacobianMatrix(config(1).JointPosition,config(2).JointPosition,config(3).JointPosition)
