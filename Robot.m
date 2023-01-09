classdef Robot
    properties
        DH_table;
        TransoformationMatrix;
        Joints;
        GeometricalJacobianMatrix;
        AnalyticalJacobianMatrix;
    end
    properties (Access = private)
        r1,inertia_vars,q, q_dot, q_dot_dot, masses, lengths, g, f_e, mu_e;
    end
    methods
        %%Default constructor
        function obj = Robot()
            obj.r1 = sym("r1"); %base frame r in z axix
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
                [0,obj.lengths(1,3),0,obj.lengths(3,3)];
                [obj.r1,0,obj.lengths(2,3)+obj.q(2),0];
                ];

             obj.f_e = [sym("fex");sym("fey");sym("fez")];
             obj.mu_e = [sym("muex");sym("muey");sym("muez")];

             obj.TransoformationMatrix = obj.getTransformationMatrixEE(obj.DH_table);

             obj.GeometricalJacobianMatrix = obj.generateJacobianMatrix(obj.DH_table,obj.Joints, 3);

             obj.AnalyticalJacobianMatrix = obj.generateAnalyticalJacobianMatrix(obj.TransoformationMatrix(1:3,4),obj.q);
        end

        function kinematics = directKinematics(obj, t_1, d_1, t_2)
           kinematics = double(subs(obj.TransoformationMatrix, [obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3),obj.q(1), obj.q(3),obj.q(2)], [0.15, 0.4, 0.3, 0.16, t_1,t_2, d_1]));
        end

        function inverKinematics = inverseKinematics(obj, px, py, pz)
            sin_theta_2 = (px/-obj.lengths(3,3));
            cos_theta_2 = sqrt(1 - power(sin_theta_2,2));
            atan_theta_2 = atan2(sin_theta_2,cos_theta_2);
            theta_2_val = double(real(subs(atan_theta_2,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3)],[0.15,0.4,0.3,0.16])));

            r = obj.lengths(1,3) + obj.lengths(3,3)*cos(theta_2_val);
            d_1 = sqrt(power(py,2) + power((pz - obj.r1),2) - power(r,2)) - obj.lengths(2,3);
            d_1_val = double(real(subs(d_1,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3)],[0.15,0.4,0.3,0.16])));

            cos_theta_1 = ((pz-obj.r1)*(obj.lengths(2,3)+d_1_val)+ r*py)/(power((obj.lengths(2,3)+d_1_val),2)+power(r,2));
            sin_theta_1 = -1*((py*(obj.lengths(2,3)+d_1_val))+(r*obj.r1)-(r*pz))/(power(r,2)+power((obj.lengths(2,3)+d_1_val),2));
            atan_theta_1 = atan2(sin_theta_1,cos_theta_1);
            theta_1_val = double(real(subs(atan_theta_1,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3)],[0.15,0.4,0.3,0.16])));
            
            inverKinematics = [theta_1_val*(180/pi), d_1_val, theta_2_val*(180/pi)];
        end

        function JacobianMatrixSolved = solveGeometricalJacobianMatrix(obj, t_1,d_1,t_2)
            JacobianMatrixSolved = double(subs(obj.GeometricalJacobianMatrix,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3), obj.q(1),obj.q(2),obj.q(3)],[0.15,0.4,0.3,0.16,t_1,d_1,t_2]));
        end

        function AnalyticalJacobianMatrixSolved = solveAnalyticalJacobianMatrix(obj, t_1,d_1,t_2)
            AnalyticalJacobianMatrixSolved = double(subs(obj.AnalyticalJacobianMatrix,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3), obj.q(1),obj.q(2),obj.q(3)],[0.15,0.4,0.3,0.16,t_1,d_1,t_2]));
        end

        function torque = langrangianEquationOfMotion(obj)
            [KE,Bq] = obj.kineticEnergy(obj.masses, obj.GeometricalJacobianMatrix, obj.q_dot, obj.DH_table,obj.Joints,obj.lengths );
            PE = obj.potentialEnergy(obj.masses, obj.g,obj.DH_table);
            Cq = obj.coriolisMatrix(obj.q,obj.q_dot, Bq);
            Gq = obj.gravityMatrix(obj.masses, obj.g, obj.GeometricalJacobianMatrix);

            torque = Bq*obj.q_dot_dot + Cq*obj.q_dot + Gq;
            simplify(torque);
        end

        function torque = newtonEulerEquationOfMotion(obj)
            [w, w_dot, p_dot_dot_c] = obj.forwardRecursion(obj.q, obj.q_dot, obj.q_dot_dot, obj.g, obj.Joints, obj.lengths);
            torque = obj.backwardEquation(obj.masses,p_dot_dot_c, w, w_dot, obj.f_e, obj.mu_e, obj.Joints, obj.lengths );
        end

        function test(obj, robot)
           q_values = [2.357, 0.1498, 1.1815];

           density = 2700; %Almunium density
           mass_values = [obj.lengths(1,1)*obj.lengths(1,2)*obj.lengths(1,3)*density, obj.lengths(2,1)*obj.lengths(2,2)*obj.lengths(3,3)*density, obj.lengths(3,1)*obj.lengths(3,2)*obj.lengths(3,3)*density];
           lenghts_values = [[0.02, 0.1256, 0.4]; [0.03, 0.03, 0.3]; [0.02, 0.1256, 0.16]];

           q_dot_values = [0,0,0];
           q_dot_dot_values = [0,0,0];
           
           g_value = 9.806;

           f_e_values = [0;0;0];
           mu_e_values = [0;0;0];

           disp("q vector: ");
           disp(q_values);

           config = homeConfiguration(robot);
           config(1).JointPosition = q_values(1);
           config(2).JointPosition = q_values(2);
           config(3).JointPosition = q_values(3);

           transformResultToolKit = getTransform(robot, config, "ee");
           transformResultOur = obj.directKinematics(q_values(1),q_values(2),q_values(3));

           disp("Toolkit Direct Kinematics: ");
           disp(transformResultToolKit);

           disp("Ours Direct Kinematics: ");
           disp(transformResultOur);

           disp("Inverse Kinematics q's: ");
           transformationRes = obj.directKinematics(q_values(1)*(pi/180),q_values(2),q_values(3)*(pi/180));
           disp(obj.inverseKinematics(transformationRes(1,4),transformationRes(2,4),transformationRes(3,4)));

           disp("Geometrical Jacobian Toolkit: ");
           disp(geometricJacobian(robot, config,"ee"));

           disp("Geometrical Jacobian Ours: ");
           disp(obj.solveGeometricalJacobianMatrix(q_values(1),q_values(2),q_values(3)));

           disp("Analytical Jacobian Ours: ");
           disp(obj.solveAnalyticalJacobianMatrix(q_values(1),q_values(2),q_values(3)));

           disp("<--B(q) validation-->");
           [~,Bq] = obj.kineticEnergy(obj.masses, obj.GeometricalJacobianMatrix, obj.q_dot, obj.DH_table,obj.Joints,obj.lengths);
           Bq_subs = Bq;
           for i = 1:size(obj.Joints,2)
               Bq_subs = subs(Bq_subs, [obj.q(i), obj.q_dot(i), obj.masses(i)],[q_values(i),q_dot_values(i), mass_values(i)]);
               Bq_subs = subs(Bq_subs, [obj.lengths(i,1),obj.lengths(i,2),obj.lengths(i,3)],[lenghts_values(i,1),lenghts_values(i,2),lenghts_values(i,3)]);
           end

           disp("B(q) values: ");
           Bq_solved = double(Bq_subs);
           disp(Bq_solved);

           if(Bq_solved' == Bq_solved)
               disp("B(q) is skew-symetric matrix.");
           else
               disp("B(q) is not skew-symetric matrix.");
           end

           disp("Eigen vector: ");
           disp(eig(Bq_solved));

           A = Bq_solved;
           eig_vals = eig((A + A')/2); 
           is_pd = [sum([sign(eig_vals)==1])==length(A)];  % flag to check if it is PD
           if is_pd ==1 
                disp("A is Positive Definite.")
           else 
                disp("A isn't Positive Definite.")
           end

           langrangian_EOM = obj.langrangianEquationOfMotion();

           for i = 1:size(obj.Joints,2)
               langrangian_EOM = subs(langrangian_EOM, [obj.q(i), obj.q_dot(i),obj.q_dot_dot(i), obj.masses(i)],[q_values(i),q_dot_values(i),q_dot_dot_values(i), mass_values(i)]);
               langrangian_EOM = subs(langrangian_EOM, [obj.lengths(i,1),obj.lengths(i,2),obj.lengths(i,3)],[lenghts_values(i,1),lenghts_values(i,2),lenghts_values(i,3)]);
           end

           langrangian_EOM = subs(langrangian_EOM, obj.g,g_value);
           disp("Tourqes Langrangian: ");
           tourqesLangrangian = double(langrangian_EOM)

           newtonEuler_EOM = obj.newtonEulerEquationOfMotion();
           for i = 1:size(obj.Joints,2)
               newtonEuler_EOM = subs(newtonEuler_EOM, [obj.q(i), obj.q_dot(i),obj.q_dot_dot(i), obj.masses(i)],[q_values(i),q_dot_values(i),q_dot_dot_values(i), mass_values(i)]);
               newtonEuler_EOM = subs(newtonEuler_EOM, [obj.lengths(i,1),obj.lengths(i,2),obj.lengths(i,3)],[lenghts_values(i,1),lenghts_values(i,2),lenghts_values(i,3)]);
           end

           disp("Tourqes Newton Euler: ");
           newtonEuler_EOM = subs(newtonEuler_EOM, [obj.g, obj.f_e(1),obj.f_e(2),obj.f_e(3), obj.mu_e(1),obj.mu_e(2),obj.mu_e(3)],[g_value, f_e_values(1),f_e_values(2),f_e_values(3),mu_e_values(1),mu_e_values(2),mu_e_values(3)]);
           tourqesNewtonEuler = double(newtonEuler_EOM)

        end

        function demo(obj)
            
           
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

            if(index == (size(obj.DH_table, 2)))
                rotationHomogeneousMatrix = [
                [0,0,1,0];
                [0,1,0,0];
                [-1,0,0,0];
                [0,0,0,1];
            ];

            matrixA = matrixA * rotationHomogeneousMatrix;
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
            analyticalJacobianMatrix = sym(zeros(3,size(q,1)));
            
            for j = 1:3
               for i = 1:size(q,1)
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
        function inertia = momentOfInertia(~,Bq, q_dot_dot, i,n)
            inertia = 0;
            for j = 1:n
                inertia = Bq(i,j)*q_dot_dot;
            end
        end

        function Coriolis = coriolisMatrix(~, q, q_dot, Bq)
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
        function Gravity = gravityMatrix(~, m, g, J)
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

        function [w_app,w_dot_app,p_dot_dot_c_app] = forwardRecursion(obj,q,q_dot,q_dot_dot, g, joints, lengths)
            w_pre = [0;0;0];
            w_pre_dot = [0;0;0];
            Zo = [0;0;1];
            w = 0;
            w_dot = 0;
            p_dot_dot_pre = [0;0;-g];

            p_dot_dot_c_app = sym(zeros(3,1,size(q,1)));
            w_app = sym(zeros(3,1,size(q,1)));
            w_dot_app = sym(zeros(3,1,size(q,1)));

            r_c = [
                [-lengths(1,3)/2,0,0];
                [0,0,-lengths(2,3)/2];
                [0,0,-lengths(3,3)/2];
                ];

            for i = 1:size(joints,2)
                H = (obj.dhMatrixGenerator(obj.DH_table, i)*obj.dhMatrixGenerator(obj.DH_table, i+1));
                R = H(1:3,1:3);
                
                w = R'*w_pre;
                w_dot = R'*w_pre_dot;

                if(joints(i) == "Revolut")
                    w = w + (R'*q_dot(i)*Zo);
                    w_dot = w_dot + (R'*(q_dot_dot(i)*Zo + q_dot(i)*cross(w_pre,Zo)));
                end
                r_l = R'*H(1:3,4); %link i-1 to i
                p_dot_dot = R'*p_dot_dot_pre+cross(w_dot,r_l)+cross(w,cross(w, r_l));

                if(joints(i) == "Prismatic")
                    p_dot_dot = p_dot_dot + (R'*q_dot_dot(i)*Zo + cross( 2*q_dot(i)*w, R'*Zo));
                end

                p_dot_dot_c = p_dot_dot + cross(w_dot, r_c(i,1:3)') + cross(w, cross(w,r_c(i,1:3)'));

                p_dot_dot_c_app(:,:,i) = p_dot_dot_c;
                w_app(:,:,i) = w;
                w_dot_app(:,:,i) = w_dot;

                w_pre = w;
                p_dot_dot_pre = p_dot_dot;
            end
        end

        function tourqes = backwardEquation(obj, masses, p_dot_dot_c,w,w_dot, f_e, mu_e, joints, lengths)
            he = [f_e;mu_e];
            f_next = f_e;
            mu_next = mu_e;
            tourqes = sym(zeros(size(joints,2),1));
            Zo = [0;0;1];

            r_c = [
                [-lengths(1,3)/2,0,0];
                [0,0,-lengths(2,3)/2];
                [0,0,-lengths(3,3)/2];
                ];

            for i = size(joints,2):-1:1
                H = (obj.dhMatrixGenerator(obj.DH_table, i)*obj.dhMatrixGenerator(obj.DH_table, i+1));
                R = H(1:3,1:3);

                r_l = R'*H(1:3,4); %link i-1 to i
                
                f = R*f_next+masses(i)*p_dot_dot_c(:,:,i);

                if joints(i) == "Revolut"
                    Isym = obj.cylinderInertia(lengths(i,1), lengths(i,2), lengths(i,3), masses(i));
                elseif joints(i) == "Prismatic"
                    Isym = obj.rectangleInertia(lengths(i,1),lengths(i,2), lengths(i,3), masses(i));
                end
                I = obj.translateInertia(Isym, masses(i), (-lengths(i,3)/2));
                mu = cross(-f, r_l+r_c(i,1:3)') + R*mu_next + cross(R*f_next, r_c(i,1:3)') + I*w_dot(:,:,i) + cross(w(:,:,i),I*w(:,:,i));

                 if joints(i) == "Revolut"
                    tourqes(i) = f'*R'*Zo;
                elseif joints(i) == "Prismatic"
                    tourqes(i) = mu'*R'*Zo;
                 end

                 mu_next = mu;
                 f_next = f;
            end
        end
    end
end