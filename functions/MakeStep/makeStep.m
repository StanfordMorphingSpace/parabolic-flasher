function [p_new, v_new, E_v] = makeStep(p, v, a, getJacobian, getb, labels, beta, i, mass, dt, varargin)
    a_stack = a'; a_stack = a_stack(:);
    v_stack = v'; v_stack = v_stack(:);
    p_stack = p'; p_stack = p_stack(:);

    J = getJacobian(p_stack, labels, beta, varargin);

    J_pseudo = J'/(J*J');

    proj = (speye(length(a_stack)) - J_pseudo*J);

    r_stack = a_stack'*proj; r_stack = r_stack';

    if i ~= 1
        v_stack = proj*v_stack.*norm(v_stack)/norm(proj*v_stack);
    end

    v_stack = v_stack + r_stack.*dt;
    p_stack = p_stack + v_stack.*dt;

    % pullback
    b = getb(p_stack, labels, beta, varargin);
    for j = 1:100
        if norm(b)<1e-6
            break
        end
        dp_stack = -J_pseudo*b;
        p_stack = p_stack + dp_stack;

        J = getJacobian(p_stack, labels, beta, varargin);
        J_pseudo = J'/(J*J');

        b = getb(p_stack, labels, beta, varargin);
    end

    p_new = reshape(p_stack, size(p')); p_new = p_new';
    v_new = reshape(v_stack, size(v')); v_new = v_new';
    E_v = 1/2.*mass.*vecnorm(v')'.^2;
end
