function w_optimal = eigen_decomposition(W_optimal, N, H_E, H_ B1, constraint, P)
    [eigenVectors_matrix, eigenValues_matrix] = eig(W_optimal);
    
    [~, max_index] = max(diag(eigenValues_matrix));
    w_l = sqrt(eigenValues_matrix(max_index, max_index)) * eigenVectors_matrix(:, max_index);
    
    % Check whether w_optimal satisfies Constraint 1
    tmp1 = norm(H_B1 * w_l)^2 - constraint;
    
    if tmp1 < 0
        scale = sqrt(constraint) / norm(H_B1 * w_l) + 0.01;
        w_scale = w_l * scale;
    else
        w_scale = w_l;
    end
    
    % Check Constraint 1 violation
    if norm(H_B1 * w_scale)^2 - constraint < 0
        flag1 = 1; % Constraint 1 is violated
    end
    
    % Check Constraint 2 violation
    if norm(w_scale)^2 - P > 0
        flag2 = 1; % Constraint 2 is violated
    end

    w_optimal = w_scale;
end