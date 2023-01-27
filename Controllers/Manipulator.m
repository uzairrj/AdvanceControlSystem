function [sys, x, str, ts] = Manipulator(t,x,u,flag)
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
        B = double(inertiaMatrix(q));
        C = double(corolisMatrix(q,q_dot));
        G = double(gravityMatrix(q));
        
        q_dot_dot = inv(B)*(tau - C*q_dot - G);

        sys = [q_dot_dot; q_dot];

    elseif flag == 3 %generate the output
        sys = x;
    end
end

function B = inertiaMatrix(q)
    q1 = q(1); % not used in this matrix
    q2 = q(2);
    q3 = q(3);

    B = [[(1377*q2)/6250 + (5286372695289*pi)/76293945312500 + (864*cos(q3))/78125 + (864*cos(q3)^2)/390625 + (351*q2^2)/625 + 782036240702373/6103515625000000, (216*cos(q3))/15625 + 702/3125, (216*sin(q3)*(q2 + 3/10))/15625]
        [                                                                                                                       (216*cos(q3))/15625 + 702/3125,                        351/625,                               0]
        [                                                                                                                      (216*sin(q3)*(q2 + 3/10))/15625,                              0,                881199/244140625]];
end

function C = corolisMatrix(q, q_dot)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);

    q1_dot = q_dot(1);
    q2_dot = q_dot(2);
    q3_dot = q_dot(3);

    C = [[q2_dot*((351*q2)/625 + 1377/12500) - q3_dot*((432*sin(2*q3))/390625 + (432*sin(q3))/78125), q1_dot*((351*q2)/625 + 1377/12500), (216*q3_dot*cos(q3)*(q2 + 3/10))/15625 - q1_dot*((432*sin(q3))/78125 + (864*cos(q3)*sin(q3))/390625)]
        [                         -q1_dot*((351*q2)/625 + 1377/12500) - (216*q3_dot*sin(q3))/15625,                                  0,                                                                            -(216*q1_dot*sin(q3))/15625]
        [        q1_dot*((432*sin(2*q3))/390625 + (432*sin(q3))/78125) + (216*q2_dot*sin(q3))/15625,         (216*q1_dot*sin(q3))/15625,                                                                                                    0]];
end

function G = gravityMatrix(q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    G = [(1720953*cos(q1))/781250 - (6751431*sin(q1))/6250000 + (264762*cos(q1)*cos(q3))/1953125 + (132381*pi*cos(q1))/156250 - (1720953*q2*sin(q1))/312500;
                                                                                                                          (1720953*cos(q1))/312500;
                                                                                                                 -(264762*sin(q1)*sin(q3))/1953125];
end

