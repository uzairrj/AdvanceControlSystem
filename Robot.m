classdef Robot
    properties
        DH_table;
        TransoformationMatrix;
        Joints;
        GeometricalJacobianMatrix;
        AnalyticalJacobianMatrix;
    end
    properties (Access = private)
        r1,inertia_vars,q, q_dot, q_dot_dot, masses, lengths, g, f_e, mu_e, P_COM;
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
            obj.P_COM = [
                [-obj.lengths(1,3)/2,0,0];
                [0,0,-obj.lengths(2,3)/2];
                [0,0,-obj.lengths(3,3)/2];]

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

             obj.TransoformationMatrix = obj.getTransformationMatrix(obj.DH_table, 4);

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

        function torque = newtonEulerEquationOfMotion(obj)
            [w, w_dot, p_dot_dot_c] = obj.forwardRecursion(obj.q, obj.q_dot, obj.q_dot_dot, obj.g, obj.Joints, obj.lengths);
            torque = obj.backwardEquation(obj.masses,p_dot_dot_c, w, w_dot, obj.f_e, obj.mu_e, obj.Joints, obj.lengths );
        end

        function test(obj, robot)
           q_values = [2.357, 0.1498, 1.1815];
           %q_values = [0, 0, 0];
           r1_value = 0.15;
           density = 2700; %Almunium density
           mass_values = [obj.lengths(1,1)*obj.lengths(1,2)*obj.lengths(1,3)*density, obj.lengths(2,1)*obj.lengths(2,2)*obj.lengths(3,3)*density, obj.lengths(3,1)*obj.lengths(3,2)*obj.lengths(3,3)*density];
           
           lenghts_values = [
               [0.02, 0.1256, 0.4]; 
               [0.03, 0.03, 0.3]; 
               [0.02, 0.1256, 0.16]
               ];

           mass_values = subs(mass_values, [obj.lengths(1,1),obj.lengths(1,2), obj.lengths(1,3), obj.lengths(2,1),obj.lengths(2,2), obj.lengths(2,3),obj.lengths(3,1),obj.lengths(3,2), obj.lengths(3,3)], [lenghts_values(1,1), lenghts_values(1,2), lenghts_values(1,3),lenghts_values(2,1), lenghts_values(2,2), lenghts_values(2,3),lenghts_values(3,1), lenghts_values(3,2), lenghts_values(3,3)]);

           q_dot_values = [0,0,0];
           q_dot_dot_values = [0,0,0];
           
           g_value = 9.806;

           f_e_values = [0;0;0];
           mu_e_values = [0;0;0];

           parameters_old = [obj.masses(1), obj.masses(2), obj.masses(3),obj.g, obj.q(1), obj.q(2), obj.q(3),obj.r1, obj.lengths(1,1),obj.lengths(1,2), obj.lengths(1,3), obj.lengths(2,1),obj.lengths(2,2), obj.lengths(2,3),obj.lengths(3,1),obj.lengths(3,2), obj.lengths(3,3), obj.q_dot(1),obj.q_dot(2),obj.q_dot(3)];
           parameters_new  = [mass_values(1), mass_values(2), mass_values(3), g_value, q_values(1),q_values(2),q_values(3),r1_value, lenghts_values(1,1), lenghts_values(1,2), lenghts_values(1,3),lenghts_values(2,1), lenghts_values(2,2), lenghts_values(2,3),lenghts_values(3,1), lenghts_values(3,2), lenghts_values(3,3),q_dot_values(1),q_dot_values(2),q_dot_values(3)];

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

           disp("<<-- Testing Homogenous Transformations -->>")

           H1_1_ours = obj.getHomogeneousMatrix(obj.DH_table, 1,1);
           H1_1_tool = getTransform(robot, config, "base_link");
           disp("Toolbox H_1_1 results: ");
           disp(H1_1_tool);
           disp("Ours H_1_1 results: ");
           disp(double(subs(H1_1_ours, parameters_old,parameters_new)));

           H1_2_ours = obj.getHomogeneousMatrix(obj.DH_table, 2,3);
           H1_2_handMade = obj.dhMatrixGenerator(obj.DH_table, 2);
           disp("Handmade H_2_3 results: ");
           disp(double(subs(H1_2_handMade, parameters_old,parameters_new)));
           disp("Ours H_2_3 results: ");
           disp(double(subs(H1_2_ours, parameters_old,parameters_new)));

           H1_3_ours = obj.getHomogeneousMatrix(obj.DH_table, 1,3);
           H1_3_handMade = obj.dhMatrixGenerator(obj.DH_table, 1) * obj.dhMatrixGenerator(obj.DH_table, 2);
           disp("Handmade H_1_3 results: ");
           disp(double(subs(H1_3_handMade, parameters_old,parameters_new)));
           disp("Ours H_1_3 results: ");
           disp(double(subs(H1_3_ours, parameters_old,parameters_new)));

           H4_5_ours = obj.getHomogeneousMatrix(obj.DH_table, 4,5);
           H4_5_handMade = obj.dhMatrixGenerator(obj.DH_table, 4) * [
                [0,0,1,0];
                [0,1,0,0];
                [-1,0,0,0];
                [0,0,0,1];
            ];
           disp("Handmade H_4_5 results: ");
           disp(double(subs(H4_5_handMade, parameters_old,parameters_new)));
           disp("Ours H_4_5 results: ");
           disp(double(subs(H4_5_ours, parameters_old,parameters_new)));

           H_toolbox_ee = getTransform(robot, config, "ee");
           H_ours_ee = obj.getHomogeneousMatrix(obj.DH_table,1,5);

           disp("Toolbox H_1_5 (ee) results: ");
           disp(double(subs(H_toolbox_ee, parameters_old,parameters_new)));
           disp("Ours H_1_5(ee) results: ");
           disp(double(subs(H_ours_ee, parameters_old,parameters_new)));

           disp("<<-- Testing Partials Jacobians -->>")
           disp("J_P_2_3 Handmade: ");
           H = obj.getHomogeneousMatrix(obj.DH_table, 2,3);
           P_C_2_3 = H(1:3,1:3) * obj.P_COM(1,:)' + H(1:3,4);
           
           H_pre = obj.getHomogeneousMatrix(obj.DH_table,2,2); %H_j-1
           P_J_m_1 = H_pre(1:3,4); %P_j-1
           Z2 = H_pre(1:3,3);
           partial = sym(zeros(6,3));
           partial(4:6,1) = cross(Z2, P_C_2_3 - P_J_m_1);
           partial(1:3,1) = Z2;
           disp(simplify(partial));

           disp("J_P_2_3 Ours");
           partial_ours = obj.partialJacobian(obj.DH_table,obj.P_COM,obj.Joints, 2);
           disp(simplify(partial_ours));

           disp("J_P_2_4 Handmade: ");
           partial = sym(zeros(6,3));

           H = obj.getHomogeneousMatrix(obj.DH_table, 2,4);
           P_C_2_4 = H(1:3,1:3) * obj.P_COM(2,:)' + H(1:3,4);
           
           H_2_2 = obj.getHomogeneousMatrix(obj.DH_table,2,2); %H_j-1
           P_2_2 = H_2_2(1:3,4); %P_j-1
           Z2 = H_pre(1:3,3);

           
           partial(4:6,1) = cross(Z2, P_C_2_4 - P_2_2);
           partial(1:3,1) = Z2;

           H_2_3 = obj.getHomogeneousMatrix(obj.DH_table,2,3); %H_j-1
           partial(4:6,2) = H_2_3(1:3,3);

           disp(simplify(partial));

           disp("J_P_2_4 Ours");
           partial_ours = obj.partialJacobian(obj.DH_table,obj.P_COM,obj.Joints, 3);
           disp(simplify(partial_ours));


           disp("J_P_2_5 Handmade: ");
           partial = sym(zeros(6,3));

           H = obj.getHomogeneousMatrix(obj.DH_table, 2,5);
           P_C_2_5 = H(1:3,1:3) * obj.P_COM(3,:)' + H(1:3,4);
           
           H_2_2 = obj.getHomogeneousMatrix(obj.DH_table,2,2); %H_j-1
           P_2_2 = H_2_2(1:3,4); %P_j-1
           Z2 = H_pre(1:3,3);

           
           partial(4:6,1) = cross(Z2, P_C_2_5 - P_2_2);
           partial(1:3,1) = Z2;

           H_2_3 = obj.getHomogeneousMatrix(obj.DH_table,2,3); %H_j-1
           partial(4:6,2) = H_2_3(1:3,3);

           H_2_4 = obj.getHomogeneousMatrix(obj.DH_table,2,4); %H_j-1
           P_2_4 = H_2_4(1:3,4); %P_j-1
           partial(4:6,3) = cross(H_2_4(1:3,3), P_C_2_5 - P_2_4);
           partial(1:3,3) = H_2_4(1:3,3);

           disp(simplify(partial));

           disp("J_P_2_5 Ours");
           partial_ours = obj.partialJacobian(obj.DH_table,obj.P_COM,obj.Joints, 4);
           disp(simplify(partial_ours));

           disp("<<-- Potential Energy -->>");
           PE = obj.potentialEnergy(obj.masses, obj.DH_table,obj.P_COM ,obj.g);
           disp(double(subs(PE, parameters_old, parameters_new)));

           disp(" <<-- Gravity Matrix -->>");
           disp("Gravity Differential Matrix");
           gravityDifferential = obj.gravityMatrixDrifferential(PE, obj.q);
           disp(double(subs(gravityDifferential, parameters_old, parameters_new)));

           disp("Gravity Matrix");
           gravityPotential = obj.gravityMatrix(obj.DH_table, obj.P_COM, obj.Joints, obj.masses, obj.g);
           disp(double(subs(gravityPotential,  parameters_old, parameters_new)));

           disp("<<-- Inertia Matrix -->>");
           Bq = obj.inertiaMatrix(obj.DH_table, obj.masses, obj.P_COM, obj.Joints, obj.lengths);
           Bq_solve = double(subs(Bq, parameters_old, parameters_new));
           disp(Bq_solve);

           if(Bq_solve == Bq_solve')
               disp("Bq matrix is skew-symetric!");
           else
               disp("Bq matrix is not skew-symetric!");
           end

           eig_vals = eig(Bq_solve);
           disp("Eigen Values:");
           disp(eig_vals);

           is_pd = [sum([sign(eig_vals)==1])==length(Bq_solve)];
           if is_pd ==1 
                disp("Bq is Positive Definite")
           else 
                disp("Bq isn't Positive Definite")
           end

           disp(" <<-- Kinetic Energy -->>");
           KE = obj.kineticEnergy(Bq, obj.q_dot);
           disp(double(subs(KE, parameters_old, parameters_new)));
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

        %Matrix from 1 to index, 1 is base frame
        function matrix = getTransformationMatrix(obj, dh_table, index)
            matrix = obj.dhMatrixGenerator(dh_table,1);
            for i = 2:index
                matrix = matrix * obj.dhMatrixGenerator(dh_table,i);
            end

            if(index == (size(obj.DH_table, 2)))
                rotationHomogeneousMatrix = [
                [0,0,1,0];
                [0,1,0,0];
                [-1,0,0,0];
                [0,0,0,1];
            ];

            matrix = matrix * rotationHomogeneousMatrix;
            end
        end

        function jacobianMatrix = generateJacobianMatrix(obj, dh_table,joints, index)
            jacobianMatrix = sym(zeros(6,size(joints,2)));
            pn = obj.getTransformationMatrix(dh_table,index+1);
            pn =  pn(1:3,4); %extracting only position part
            matrix = obj.dhMatrixGenerator(dh_table,1);
            for i = 1:index %1 is base frame
                if joints(i) == "Revolut"
                    jacobianMatrix(1:3,i) = matrix(1:3,3);
                    jacobianMatrix(4:6,i) = cross(matrix(1:3,3),pn - matrix(1:3,4));
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

        %1 is the base frame
        function H = getHomogeneousMatrix(obj, dh_table, startIndex, endIndex)
            H = sym(eye(4,4));

            if(startIndex == endIndex) %if no transformation is needed, use empty homogeneous transformation
                return
            end

            for i = startIndex: endIndex-1
                H = H * obj.dhMatrixGenerator(dh_table, i);
            end

            if(endIndex == 5) %one extra rotation in the end
                 rotationHomogeneousMatrix = [
                [0,0,1,0];
                [0,1,0,0];
                [-1,0,0,0];
                [0,0,0,1];
            ];

            H = H * rotationHomogeneousMatrix;
            end
        end

        %index 1 is base, but in this the base is transfered to link 1 so
        %it is 2
        function partial = partialJacobian(obj, dh_table, P_COM, joints, index)
            partial = sym(zeros(6, size(joints,2)));
            H = obj.getHomogeneousMatrix(dh_table, 2, index+1);
            P_li = H(1:3,1:3)*P_COM(index-1,:)'+H(1:3,4);

            for i=2:index
                H_pre = obj.getHomogeneousMatrix(dh_table,2, i);

                if(joints(i-1) == "Revolut")
                    partial(4:6,i-1) = cross(H_pre(1:3,3), P_li - H_pre(1:3,4));
                    partial(1:3,i-1) = H_pre(1:3,3);
                else
                    partial(4:6,i-1) = H_pre(1:3,3);
                end

            end
        end
    
        function PE = potentialEnergy(obj, masses, dh_table, CoM, g)
            PE = 0;
            gvec = [0;-g;0];
            for i = 2:4
                H = obj.getHomogeneousMatrix(dh_table, 2,i+1);
                Pl = H(1:3,1:3)*CoM(i-1,:)'+H(1:3,4);
                PE = PE + (-masses(i-1)*gvec'*Pl);
            end
        end

        function G = gravityMatrixDrifferential(obj, U, q)
            G = [diff(U,q(1)); diff(U,q(2)); diff(U,q(3))];
        end

        function G = gravityMatrix(obj, dh_table, P_COM, joints, masses, g)
            G = sym(zeros(3,1));
            gvec = [0;-g;0];
            for i = 2:4
                for j = 2:4
                    Jp = obj.partialJacobian(dh_table,P_COM, joints, j);
                    G(i-1,1) = G(i-1,1) + (-masses(j-1)*gvec'*Jp(4:6, i-1));
                end
            end
        end
    
        function Bq = inertiaMatrix(obj, dh_table, masses, P_COM, joints, lengths)
            Bq = 0;
            for i=2:4
                H = obj.getHomogeneousMatrix(dh_table, 2, i+1);
                J = obj.partialJacobian(dh_table, P_COM, joints, i);
                if(joints(i-1) == "Revolut")
                    Ic = obj.inertiaCylender(lengths(i-1,1),lengths(i-1,2),lengths(i-1,3),masses(i-1));
                else
                    Ic = obj.inertiaCube(lengths(i-1,1),lengths(i-1,2),lengths(i-1,3),masses(i-1));
                end
                I = obj.translateInertia(Ic, masses(i-1),P_COM(i-1,:)');
                Bq = Bq + (masses(i-1)*J(4:6,:)'*J(4:6,:)+ J(1:3,:)'*H(1:3,1:3)*I*H(1:3,1:3)'*J(1:3,:));
            end
        end

        function Ic = inertiaCylender(obj, a,b,c,m)
            Ic = [
                [1/2*m*(a^2+b^2),0,0];
                [0,1/2*m*(3*(a^2+b^2)^2 + c^2),0];
                [0,0,1/2*m*(3*(a^2+b^2)^2 + c^2)]
            ];
        end

        function Ic = inertiaCube(obj, a,b,c,m)
            Ic = [
                [1/12*m*(b^2+a^2),0,0];
                [0,1/12*m*(a^2+c^2),0];
                [0,0,1/12*m*(c^2+b^2)]
            ];
        end

        function I = translateInertia(obj, Ic, m,r)
            I = Ic + m*(r'*r*eye(3,3)-r*r');
        end
        
        function KE = kineticEnergy(obj, Bq, q_dot)
            KE = q_dot'*Bq*q_dot;
        end
    end
end