classdef Robot
    properties
        DH_table;
        TransoformationMatrix;
        Joints;
        GeometricalJacobianMatrix;
        AnalyticalJacobianMatrix;
    end
    properties (Access = private)
        a1,a2,a3,a4,t_1,t_2,d_1;
    end
    methods
        %%Default constructor
        function obj = Robot()
            obj.a1 = sym("a1");
            obj.a2 = sym("a2");
            obj.a3 = sym("a3");
            obj.a4 = sym("a4");
            obj.t_1 = sym("t_1");
            obj.t_2 = sym("t_2");
            obj.d_1 = sym("d_1");

            obj.Joints = ["Revolut","Prismatic","Revolut"];

            obj.DH_table = [
                ["theta", "alpha","r","d"]
                [pi/2,obj.t_1, 0, obj.t_2];
                [pi/2,-pi/2,0,0];
                [0,obj.a2,0,obj.a4];
                [obj.a1,0,obj.a3+obj.d_1,0];
                ];

             obj.TransoformationMatrix = obj.getTransformationMatrix(obj.DH_table);

             obj.GeometricalJacobianMatrix = obj.generateJacobianMatrix(obj.DH_table,obj.Joints);

             obj.AnalyticalJacobianMatrix = obj.generateAnalyticalJacobianMatrix(obj.TransoformationMatrix(1:3,4),[obj.t_1,obj.d_1,obj.t_2]);
        end

        function kinematics = directKinematics(obj, t_1, d_1, t_2)
           kinematics = double(subs(obj.TransoformationMatrix, [obj.a1,obj.a2,obj.a3,obj.a4,obj.t_1, obj.t_2,obj.d_1], [0.15, 0.4, 0.3, 0.16, t_1,t_2, d_1]));
        end

        function inverKinematics = inverseKinematics(obj, px, py, pz)
            sin_theta_2 = (px/-obj.a4);
            cos_theta_2 = sqrt(1 - power(sin_theta_2,2));
            atan_theta_2 = atan2(sin_theta_2,cos_theta_2);
            theta_2_val = double(real(subs(atan_theta_2,[obj.a1,obj.a2,obj.a3,obj.a4],[0.15,0.4,0.3,0.16])));

            r = obj.a2 + obj.a4*cos(theta_2_val);
            d_1 = sqrt(power(py,2) + power((pz - obj.a1),2) - power(r,2)) - obj.a3;
            d_1_val = double(real(subs(d_1,[obj.a1,obj.a2,obj.a3,obj.a4],[0.15,0.4,0.3,0.16])));

            cos_theta_1 = ((pz-obj.a1)*(obj.a3+d_1_val)+ r*py)/(power((obj.a3+d_1_val),2)+power(r,2));
            sin_theta_1 = -1*((py*(obj.a3+d_1_val))+(r*obj.a1)-(r*pz))/(power(r,2)+power((obj.a3+d_1_val),2));
            atan_theta_1 = atan2(sin_theta_1,cos_theta_1);
            theta_1_val = double(real(subs(atan_theta_1,[obj.a1,obj.a2,obj.a3,obj.a4],[0.15,0.4,0.3,0.16])));
            
            inverKinematics = [theta_1_val*(180/pi), d_1_val, theta_2_val*(180/pi)];
        end

        function JacobianMatrixSolved = solveGeometricalJacobianMatrix(obj, t_1,d_1,t_2)
            JacobianMatrixSolved = double(subs(obj.GeometricalJacobianMatrix,[obj.a1,obj.a2,obj.a3,obj.a4, obj.t_1,obj.d_1,obj.t_2],[0.15,0.4,0.3,0.16,t_1,d_1,t_2]));
        end

        function AnalyticalJacobianMatrixSolved = solveAnalyticalJacobianMatrix(obj, t_1,d_1,t_2)
            AnalyticalJacobianMatrixSolved = double(subs(obj.AnalyticalJacobianMatrix,[obj.a1,obj.a2,obj.a3,obj.a4, obj.t_1,obj.d_1,obj.t_2],[0.15,0.4,0.3,0.16,t_1,d_1,t_2]));
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
        function matrixA = getTransformationMatrix(obj, dh_table)
            matrixA = obj.dhMatrixGenerator(dh_table,1);
            for i = 2:(size(dh_table, 2))
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

        function jacobianMatrix = generateJacobianMatrix(obj, dh_table,joints)
            jacobianMatrix = sym(zeros(6,size(joints,2)));
            pn = obj.TransoformationMatrix(1:3,4);
            matrix = obj.dhMatrixGenerator(dh_table,1);
            for i = 1:(size(joints,2))
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

        function analyticalJacobianMatrix = generateAnalyticalJacobianMatrix(obj,pe,q)
            analyticalJacobianMatrix = sym(zeros(3,size(q,2)));
            
            for j = 1:3
               for i = 1:size(q,2)
                   analyticalJacobianMatrix(j,i) = diff(pe(j),q(i));
               end
            end
        end
    end
end