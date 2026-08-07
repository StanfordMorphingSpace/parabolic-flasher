function [p_new, v_new, E_v] = makeStep(p, v, a, getJacobian, getb, labels, beta, i, mass, dt, varargin)
    a_stack = a'; a_stack = a_stack(:);
    v_stack = v'; v_stack = v_stack(:);
    p_stack = p'; p_stack = p_stack(:);

    J = getJacobian(p_stack, labels, beta, varargin);
    
    % project a_stack onto the nullspace of J
    x = (J*J') \ (J*a_stack);
    r_stack = a_stack - J' * x;

    if i ~= 1
        x = (J*J') \ (J*v_stack);
        Pv = v_stack - J' * x;
        v_stack = Pv * (norm(v_stack) / norm(Pv));

        % damping calc from geodesic paper
        % theta = dot(v_stack, r_stack)./norm(v_stack)./norm(r_stack);
        % gamma = 0.98*theta;
        %gamma = 0.95+theta./20;
        %gamma = theta>0;
    % else
    %     gamma = zeros(size(v_stack));
    end
    
    v_stack = v_stack + r_stack.*dt;
    %v_stack = gamma.*v_stack + r_stack.*dt;
    p_stack = p_stack + v_stack.*dt;

    % pullback
    b = getb(p_stack, labels, beta, varargin);
    for j = 1:100
        if norm(b)<1e-6
            break
        end
        x = (J*J') \ b; 

        dp_stack = -J' * x;
        p_stack = p_stack + dp_stack;

        J = getJacobian(p_stack, labels, beta, varargin);

        b = getb(p_stack, labels, beta, varargin);
    end

    p_new = reshape(p_stack, size(p')); p_new = p_new';
    v_new = reshape(v_stack, size(v')); v_new = v_new';
    E_v = 1/2.*mass.*vecnorm(v')'.^2;
end
