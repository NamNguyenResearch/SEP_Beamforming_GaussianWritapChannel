function w_optimal = rand_GaussianVector(W_optimal, N, H_E, H_B, constraint, P)
    L = 1000; % Number of randomization samples

    w_scale = zeros(N, L);
    ObjectiveValue = zeros(1, L);

    for l = 1:L
        w_l = 1/sqrt(2) * (mvnrnd(zeros(N, 1), W_optimal, 1) + 1i * mvnrnd(zeros(N, 1), W_optimal, 1));

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