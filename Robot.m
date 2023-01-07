classdef Robot
    properties
        DH_table;
        TransoformationMatrix;
        Joints;
        GeometricalJacobianMatrix;
        AnalyticalJacobianMatrix;
    end
    properties (Access = private)
        r1,r2,r3,r4, inertia_vars,q, q_dot, q_dot_dot, masses, lengths, g;
    end
    methods
        %%Default constructor
        function obj = Robot()
            obj.r1 = sym("r1");
            obj.r2 = sym("r2");
            obj.r3 = sym("r3");
            obj.r4 = sym("r4");
            obj.inertia_vars = [sym("a"), sym("b"), sym("c")];
            obj.g = sym("g");
            obj.q = [[sym("q1")]; [sym("q2")]; [sym("q3")];];
            obj.q_dot = [[sym("q1_dot")]; [sym("q2_dot")]; [sym("q3_dot")];];
            obj.q_dot_dot = [[sym("q1_dot_dot")]; [sym("q2_dot_dot")]; [sym("q3_dot_dot")];];
            obj.masses = [sym("m1"), sym("m2"), sym("m3")];
            obj.lengths = [[sym("a1"),sym("b1"),sym("c1")];
                            [sym("a2"),sym("b2"),sym("c2")];
                            [sym("a3"),sym("b3"),sym("c3")]];

            obj.Joints = ["Revolut","Prismatic","Revolut"];

            obj.DH_table = [
                ["theta", "alpha","r","d"]
                [pi/2,obj.q(1), 0, obj.q(3)];
                [pi/2,-pi/2,0,0];
                [0,obj.r2,0,obj.r4];
                [obj.r1,0,obj.r3+obj.q(2),0];
                ];

             obj.TransoformationMatrix = obj.getTransformationMatrixEE(obj.DH_table);

             obj.GeometricalJacobianMatrix = obj.generateJacobianMatrix(obj.DH_table,obj.Joints, 3);

             obj.AnalyticalJacobianMatrix = obj.generateAnalyticalJacobianMatrix(obj.TransoformationMatrix(1:3,4),obj.q);
        end

        function kinematics = directKinematics(obj, t_1, d_1, t_2)
           kinematics = double(subs(obj.TransoformationMatrix, [obj.r1,obj.r2,obj.r3,obj.r4,obj.q(1), obj.q(3),obj.q(2)], [0.15, 0.4, 0.3, 0.16, t_1,t_2, d_1]));
        end

        function inverKinematics = inverseKinematics(obj, px, py, pz)
            sin_theta_2 = (px/-obj.r4);
            cos_theta_2 = sqrt(1 - power(sin_theta_2,2));
            atan_theta_2 = atan2(sin_theta_2,cos_theta_2);
            theta_2_val = double(real(subs(atan_theta_2,[obj.r1,obj.r2,obj.r3,obj.r4],[0.15,0.4,0.3,0.16])));

            r = obj.a2 + obj.a4*cos(theta_2_val);
            d_1 = sqrt(power(py,2) + power((pz - obj.r1),2) - power(r,2)) - obj.r3;
            d_1_val = double(real(subs(d_1,[obj.r1,obj.r2,obj.r3,obj.r4],[0.15,0.4,0.3,0.16])));

            cos_theta_1 = ((pz-obj.r1)*(obj.r3+d_1_val)+ r*py)/(power((obj.r3+d_1_val),2)+power(r,2));
            sin_theta_1 = -1*((py*(obj.a3+d_1_val))+(r*obj.r1)-(r*pz))/(power(r,2)+power((obj.r3+d_1_val),2));
            atan_theta_1 = atan2(sin_theta_1,cos_theta_1);
            theta_1_val = double(real(subs(atan_theta_1,[obj.r1,obj.r2,obj.r3,obj.r4],[0.15,0.4,0.3,0.16])));
            
            inverKinematics = [theta_1_val*(180/pi), d_1_val, theta_2_val*(180/pi)];
        end

        function JacobianMatrixSolved = solveGeometricalJacobianMatrix(obj, t_1,d_1,t_2)
            JacobianMatrixSolved = double(subs(obj.GeometricalJacobianMatrix,[obj.a1,obj.a2,obj.a3,obj.a4, obj.q(1),obj.q(2),obj.q(3)],[0.15,0.4,0.3,0.16,t_1,d_1,t_2]));
        end

        function AnalyticalJacobianMatrixSolved = solveAnalyticalJacobianMatrix(obj, t_1,d_1,t_2)
            AnalyticalJacobianMatrixSolved = double(subs(obj.AnalyticalJacobianMatrix,[obj.a1,obj.a2,obj.a3,obj.a4, obj.q(1),obj.q(2),obj.q(3)],[0.15,0.4,0.3,0.16,t_1,d_1,t_2]));
        end

        function torque = langrangianEquationOfMotion(obj)
            [KE,Bq] = obj.kineticEnergy(obj.masses, obj.GeometricalJacobianMatrix, obj.q_dot, obj.DH_table,obj.Joints,obj.lengths );
            PE = obj.potentialEnergy(obj.masses, obj.g,obj.DH_table);
            Cq = obj.coriolisMatrix(obj.q,obj.q_dot, Bq);
            Gq = obj.gravityMatrix(obj.masses, obj.g, obj.GeometricalJacobianMatrix);

            torque = Bq*obj.q_dot_dot + Cq*obj.q_dot + Gq;
            simplify(torque);
        end

        function x = test(obj, a,b,c,m)
            [KE,Bq] = obj.kineticEnergy(obj.masses, obj.GeometricalJacobianMatrix, obj.q_dot, obj.DH_table,obj.Joints,obj.lengths );
            PE = obj.potentialEnergy(obj.masses, obj.g,obj.DH_table);
            I1 = obj.momentOfInertia(Bq, obj.q_dot_dot(1), 1, 3);
            %Cq = obj.coriolisMatrix(obj.q,obj.q_dot, Bq);
            Gq = obj.gravityMatrix(obj.masses, obj.g, obj.GeometricalJacobianMatrix);
        end
    end
    methods (Access = private)

        function matrix = dhMatrixGenerator(~, dh_table, i)
            matrix = [
                [cos(dh_table(2,i)), -sin(dh_table(2,i))*cos(dh_table(3,i)), sin(dh_table(2,i))*sin(dh_table(3,i)), dh_table(4,i)*cos(dh_table(2,i))];
                [sin(dh_table(2,i)), cos(dh_table(2,i))*cos(dh_table(3,i)), -cos(dh_table(2,i))*sin(dh_table(3,i)), dh_table(4,i)*sin(dh_table(2,i))];
                [0, sin(dh_table(3,i)), cos(dh_table(3,i)), dh_table(5,i)];
                [0,0,0,1];

            ];
        end

        function matrixA = getTransformationMatrix(obj, dh_table, index)
            matrixA = obj.dhMatrixGenerator(dh_table,1);
            for i = 2:index
                matrixA = matrixA * obj.dhMatrixGenerator(dh_table,i);
            end
        end

        function matrixA = getTransformationMatrixEE(obj, dh_table)
            matrixA = obj.dhMatrixGenerator(dh_table,1);
            for i = 2:(size(obj.DH_table, 2))
                matrixA = matrixA * obj.dhMatrixGenerator(dh_table,i);
            end

            rotationHomogeneousMatrix = [
                [0,0,1,0];
                [0,1,0,0];
                [-1,0,0,0];
                [0,0,0,1];
            ];

            matrixA = matrixA * rotationHomogeneousMatrix;
        end

        function jacobianMatrix = generateJacobianMatrix(obj, dh_table,joints, index)
            jacobianMatrix = sym(zeros(6,size(joints,2)));
            pn = obj.TransoformationMatrix(1:3,4);
            matrix = obj.dhMatrixGenerator(dh_table,1);
            for i = 1:index
                if joints(i) == "Revolut"
                    jacobianMatrix(1:3,i) = matrix(1:3,3);
                    jacobianMatrix(4:6,i) = cross(matrix(1:3,3),pn-matrix(1:3,4));
                else
                    jacobianMatrix(1:3,i) = [0,0,0];
                    jacobianMatrix(4:6,i) = matrix(1:3,3);
                end
                matrix = matrix * obj.dhMatrixGenerator(dh_table,i+1);
            end
        end

        function analyticalJacobianMatrix = generateAnalyticalJacobianMatrix(~,pe,q)
            analyticalJacobianMatrix = sym(zeros(3,size(q,2)));
            
            for j = 1:3
               for i = 1:size(q,2)
                   analyticalJacobianMatrix(j,i) = diff(pe(j),q(i));
               end
            end
        end

        function inertia = cylinderInertia(~, a, b, c, m)
            inertia = [
                    [1/2*m*(a^2+b^2), 0, 0];
                    [0, 1/2*m*(3*(a^2+b^2)^2+c^2), 0];
                    [0, 0, 1/2*m*(3*(a^2+b^2)^2+c^2)]
                    ];
        end

        function inertia = rectangleInertia(~,a,b,c,m)
            inertia = [
                [1/12*m*(b^2+c^2), 0, 0];
                [0,1/12*m*(a^2+c^2),0];
                [0,0,1/12*m*(a^2+b^2)];
                ];
        end
        function PE = potentialEnergy(obj, m,g, dh_table)
            gVec = [[0];[0];[-g]];
            PE = 0;
            for i = 1: size(m,2)
                T = obj.getTransformationMatrix(dh_table,i+1);
                PE = PE + (m(i)*gVec'*T(1:3,4));
            end
            PE = -1*PE;
        end
        %note, May be last rotation is not added and that can cause issue
        %in R0_E rotation
        function [KE,Bq] = kineticEnergy(obj,m,J,q_dot, dh_table,joints, lenghts)
            Jz = sym(zeros(3,size(q_dot,1)));
            %Jz(1:3,1:1) = J(1:3,1:1);
            Isym = 0;
            Bq = 0;
            for i = 1:size(q_dot,1)
                Jp = Jz;
                Jp(1:3,1:i) = J(4:6,1:i);
                KEL = m(i)*Jp'*Jp;
                R = obj.getTransformationMatrix(dh_table, i); %may be i + 1
                R = R(1:3,1:3);
                if joints(i) == "Revolut"
                    Isym = obj.cylinderInertia(lenghts(i,1), lenghts(i,2), lenghts(i,3), m(i));
                elseif joints(i) == "Prismatic"
                    Isym = obj.rectangleInertia(lenghts(i,1),lenghts(i,2), lenghts(i,3), m(i));
                end
                I = obj.translateInertia(Isym, m(i), (-lenghts(i,3)/2));
                Jo = Jz;
                Jo(1:3,1:i) = J(1:3,1:i);
                KER = Jo'*R*I*R'*Jo;
                Bq = Bq + (KEL + KER);
            end
            KE = 1/2 * q_dot' * Bq * q_dot;
        end

        function inertia = translateInertia(~, inertiaTensor,m, r)
            r_ = [r,0,0];
            inertia = inertiaTensor+(m*(r_'*r_*eye(3)-r_*r_'));
        end

        %not used, computed directly
        function inertia = momentOfInertia(obj,Bq, q_dot_dot, i,n)
            inertia = 0;
            for j = 1:n
                inertia = Bq(i,j)*q_dot_dot;
            end
        end

        function Coriolis = coriolisMatrix(obj, q, q_dot, Bq)
            Coriolis = sym(zeros(size(q,1),size(q,1)));
            for i = 1:size(q,1)
                for j = 1:size(q,1)
                    for k = 1:size(q,1)
                        ck = 1/2*(diff(Bq(i,j), q(k)) + diff(Bq(i,k), q(j)) - diff(Bq(j,k), q(i))) * q_dot(k);
                        Coriolis(i,j) =  Coriolis(i,j) + ck;
                    end
                end
            end
        end
        
        %may be it is wrong
        function Gravity = gravityMatrix(obj, m, g, J)
            Gravity = sym(zeros(size(m,2),1));
            gVec = [[0];[0];[-g]];
            Jz = sym(zeros(3,size(m,2)));
            for i = 1:size(m,2)
                for j = 1:size(m,2)
                    Jp = Jz;
                    Jp(1:3,1:j) = J(4:6,1:j); %not sure. we have to use linear or angular part (Linear part should be considered, but currently using angular part)
                    JL = gVec' * Jp(1:3,i);
                    Gravity(i) = Gravity(i) + (m(j)*JL);
                end
            end

        end
    end
end