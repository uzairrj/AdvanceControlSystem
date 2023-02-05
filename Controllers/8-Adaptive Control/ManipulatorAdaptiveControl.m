function [sys, x, str, ts] = ManipulatorAdaptiveControl(t,x,u,flag, parameters)
    if flag == 0 %initialization
        x = zeros(2,1);
        str = [];
        ts = [0 0];

        sizes = simsizes;
        sizes.NumContStates = 2; %number of continuous states
        sizes.NumDiscStates = 0; %number of discrete states
        sizes.NumOutputs = 2; % number of outputs
        sizes.NumInputs = 1; %number of inputs
        sizes.DirFeedthrough = 0; % Does u depends on the output
        sizes.NumSampleTimes  = 1; %Number of function to run in time t

        sys = simsizes(sizes);
    elseif flag == 1 %generate the derivatives
        q_dot = x(1);
        q = x(2);
        tau = u(1);

        I = parameters.I;
        F = parameters.F;
        G = parameters.G;
        
        q_dot_dot = real(I*(tau-G*sin(q)-F*q_dot));

        sys = [q_dot_dot; q_dot];

    elseif flag == 3 %generate the output
        sys = x;
    end
end

