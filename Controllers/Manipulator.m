function [sys, x, str, ts] = Manipulator(t,x,u,flag, parameters)
    if flag == 0 %initialization
        x = zeros(6,1);
        str = [];
        ts = [0 0];

        sizes = simsizes;
        sizes.NumContStates = 6; %number of continuous states
        sizes.NumDiscStates = 0; %number of discrete states
        sizes.NumOutputs = 6; % number of outputs
        sizes.NumInputs = 3; %number of inputs
        sizes.DirFeedthrough = 0; % Does u depends on the output
        sizes.NumSampleTimes  = 1; %Number of function to run in time t

        sys = simsizes(sizes);
    elseif flag == 1 %generate the derivatives
        q_dot = x(1:3);
        q = x(4:6);
        tau = u(1:3);
        parameters.q = q;
        parameters.q_dot = q_dot;
        B = double(inertiaMatrix(parameters));
        C = double(corolisMatrix(parameters));
        G = double(gravityMatrix(parameters));
        
        q_dot_dot = pinv(B)*(tau - C*q_dot - G);

        sys = [q_dot_dot; q_dot];

    elseif flag == 3 %generate the output
        sys = x;
    end
end

function B = inertiaMatrix(parameters)
    q1 = parameters.q(1); % not used in this matrix
    q2 = parameters.q(2);
    q3 = parameters.q(3);

    m1 = parameters.m(1);
    m2 = parameters.m(2);
    m3 = parameters.m(3);

    B = [[(195791581307*m1)/1220703125000 + (8503*m2)/40000 + (321279862557*m3)/1220703125000 + (8*m3*cos(q3)^2)/625 + (3*m2*q2)/10 + (3*m3*q2)/5 + m2*q2^2 + m3*q2^2 + (8*m3*cos(q3))/125, (2*m2)/5 + (2*m3)/5 + (2*m3*cos(q3))/25, (2*m3*sin(q3)*(q2 + 3/10))/25]
        [                                                                                                                                         (2*m2)/5 + (2*m3)/5 + (2*m3*cos(q3))/25,                                 m2 + m3,                             0]
        [                                                                                                                                                   (2*m3*sin(q3)*(q2 + 3/10))/25,                                       0,            (32637*m3)/1562500]];
end

function C = corolisMatrix(parameters)
    q1 = parameters.q(1); %not used in this matrix
    q2 = parameters.q(2);
    q3 = parameters.q(3);

    q1_dot = parameters.q_dot(1);
    q2_dot = parameters.q_dot(2);
    q3_dot = parameters.q_dot(3);

    m1 = parameters.m(1); %not used in this matrix
    m2 = parameters.m(2);
    m3 = parameters.m(3);


    C = [[q2_dot*((3*m2)/20 + (3*m3)/10 + m2*q2 + m3*q2) - q3_dot*((4*m3*sin(2*q3))/625 + (4*m3*sin(q3))/125), q1_dot*((3*m2)/20 + (3*m3)/10 + m2*q2 + m3*q2), (2*m3*q3_dot*cos(q3)*(q2 + 3/10))/25 - q1_dot*((4*m3*sin(q3))/125 + (8*m3*cos(q3)*sin(q3))/625)]
        [                        -q1_dot*((3*m2)/20 + (3*m3)/10 + m2*q2 + m3*q2) - (2*m3*q3_dot*sin(q3))/25,                                              0,                                                                       -(2*m3*q1_dot*sin(q3))/25]
        [                      q1_dot*((4*m3*sin(2*q3))/625 + (4*m3*sin(q3))/125) + (2*m3*q2_dot*sin(q3))/25,                       (2*m3*q1_dot*sin(q3))/25,                                                                                               0]];
end

function G = gravityMatrix(parameters)
    q1 = parameters.q(1);
    q2 = parameters.q(2);
    q3 = parameters.q(3);
    
    g = parameters.g;

    m1 = parameters.m(1);
    m2 = parameters.m(2);
    m3 = parameters.m(3);

    G = [(g*m1*cos(q1))/5 - g*m2*((3*sin(q1))/20 - (2*cos(q1))/5 + q2*sin(q1)) + g*m3*((2*cos(q1))/5 + (2*cos(q1)*cos(q3))/25 - sin(q1)*(q2 + 3/10));
                                                                                                                        g*cos(q1)*(m2 + m3);
                                                                                                               -(2*g*m3*sin(q1)*sin(q3))/25];
end

