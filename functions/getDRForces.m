function [F_axial, E_axial, F_crease, E_crease, F_damping] = getDRForces(nodes_f, nodes, edges, adj_faces, EA, k_fold, lengths, angles, gamma, v, mass)
    % nodes only used to calculate correct norm direction (based on flat)
    F_crease = zeros(length(nodes_f), 3);
    E_crease = zeros(length(edges), 1);

    F_axial = zeros(length(nodes_f), 3);
    E_axial = zeros(length(edges), 1);

    F_damping = zeros(length(nodes), 3);
    for i = 1:length(edges)
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);

        p1 = nodes_f(p1_index, :);
        p2 = nodes_f(p2_index, :);
        p1_u = nodes(p1_index, :);
        p2_u = nodes(p2_index, :);
        
        % bending forces
        if length(adj_faces{i}) == 2
            k_crease = lengths(i)*k_fold(i);
            p3_1_index = adj_faces{i}(1);
            p3_2_index = adj_faces{i}(2);

            p3_1 = nodes_f(p3_1_index, :);
            p3_2 = nodes_f(p3_2_index, :); 
            p3_1_u = nodes(p3_1_index, :);
            p3_2_u = nodes(p3_2_index, :);

            h1 = norm(cross(p1-p3_1, p2-p3_1))/norm(p2-p1);
            h2 = norm(cross(p1-p3_2, p2-p3_2))/norm(p2-p1);

            n1_u = cross(p1_u-p3_1_u, p2_u-p3_1_u); n1_u=n1_u/norm(n1_u);
            n2_u = cross(p1_u-p3_2_u, p2_u-p3_2_u); n2_u=n2_u/norm(n2_u);

            n1 = cross(p1-p3_1, p2-p3_1); n1=n1/norm(n1);
            n2 = cross(p3_2-p1, p3_2-p2); n2=n2/norm(n2);

            if n1_u(3) < 0
                n1 = -n1;
            end

            if n2_u(3) < 0
                n2 = -n2;
            end

            angle = angleBetweenVectors3d(n1, n2); angle = real(angle);

            if dot(n1, p3_2-p3_1) < 0
                angle = -angle;
            end

            a431 = angleBetweenVectors3d(p1-p2, p3_1-p2);
            a314 = angleBetweenVectors3d(p3_1-p1, p2-p1);
            a423 = angleBetweenVectors3d(p3_2-p2, p1-p2);
            a342 = angleBetweenVectors3d(p2-p1, p3_2-p1);

            dthdp1 = -cot(a431)/(cot(a314)+cot(a431)).*n1./h1 + -cot(a423)/(cot(a342)+cot(a423)).*n2./h2;
            dthdp2 = -cot(a314)/(cot(a314)+cot(a431)).*n1./h1 + -cot(a342)/(cot(a342)+cot(a423)).*n2./h2;

            F_crease(p3_1_index, :) = F_crease(p3_1_index, :) - k_crease.*(angle-angles(i)).*n1./h1;
            F_crease(p3_2_index, :) = F_crease(p3_2_index, :) - k_crease.*(angle-angles(i)).*n2./h2;
    
            F_crease(p1_index, :) = F_crease(p1_index, :) - k_crease.*(angle-angles(i)).*dthdp1;
            F_crease(p2_index, :) = F_crease(p2_index, :) - k_crease.*(angle-angles(i)).*dthdp2;

            % E_crease(p3_1_index) = E_crease(p3_1_index) + k_crease/h1*(angle-angles(i))^2;
            % E_crease(p3_2_index) = E_crease(p3_2_index) + k_crease/h2*(angle-angles(i))^2;
            % 
            % E_crease(p1_index) = E_crease(p1_index) + k_crease*norm(dthdp1)*(angle-angles(i))^2;
            % E_crease(p2_index) = E_crease(p2_index) + k_crease*norm(dthdp2)*(angle-angles(i))^2;
            E_crease(i) = (1/2)*k_crease*(angle-angles(i))^2;
        end
        
        % Stretching forces
        dir12 = (p2-p1)/norm(p1-p2);
        norm12 = norm(p1-p2);

        F_axial(p1_index, :) = F_axial(p1_index, :) + (EA/lengths(i)).*(norm12-lengths(i)).*dir12;
        F_axial(p2_index, :) = F_axial(p2_index, :) - (EA/lengths(i)).*(norm12-lengths(i)).*dir12;
        
        % E_axial(p1_index) = E_axial(p1_index) + (EA/lengths(i)).*(norm12-lengths(i))^2;
        % E_axial(p2_index) = E_axial(p2_index) + (EA/lengths(i)).*(norm12-lengths(i))^2;
        E_axial(i) = (1/2)*(EA/lengths(i)).*(norm12-lengths(i))^2;
        
        % Damping
        c1 = 2*gamma*sqrt(EA/lengths(i)*mass(p1_index, 1));
        c2 = 2*gamma*sqrt(EA/lengths(i)*mass(p2_index, 1));

        F_damping(p1_index, :) = F_damping(p1_index, :) + c1.*(v(p2_index, :)-v(p1_index, :));
        F_damping(p2_index, :) = F_damping(p2_index, :) + c2.*(v(p1_index, :)-v(p2_index, :));
    end
end