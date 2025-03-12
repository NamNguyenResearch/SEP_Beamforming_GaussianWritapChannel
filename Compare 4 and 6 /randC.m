function w_optimal = randC(W_optimal, N, H_E, H_B, constraint, P)
    L = 1000; % Number of randomization samples
    
    [eigenVectors_matrix, eigenValues_matrix] = eig(W_optimal);
    
    w_scale = zeros(N, L);
    ObjectiveValue = zeros(1, L);

    for l = 1:L
        % Generate random complex vector v_l
        v_l = 1/sqrt(2) * (randn(N, 1) + 1i * randn(N, 1));

        % Generate random beamforming vector w_l
        w_l = eigenVectors_matrix * sqrt(eigenValues_matrix) * v_l;

        % Check whether w_optimal satisfies Constraint 1
        tmp1 = norm(H_B * w_l)^2 - constraint;

        if tmp1 < 0
            scale = sqrt(constraint) / norm(H_B * w_l) + 0.01;
            w_scale(:, l) = w_l * scale;
        else
            w_scale(:, l) = w_l;
        end

        % Check Constraint 1 violation
        if norm(H_B * w_scale(:, l))^2 - constraint < 0
            flag1 = 1; % Constraint 1 is violated
        end

        % Check Constraint 2 violation
        if norm(w_scale(:, l))^2 - P > 0
            flag2 = 1; % Constraint 2 is violated
        end

        ObjectiveValue(l) = norm(H_E * w_scale(:, l))^2;
    end

    [~, index_min] = min(ObjectiveValue);
    w_optimal = w_scale(:, index_min);
end