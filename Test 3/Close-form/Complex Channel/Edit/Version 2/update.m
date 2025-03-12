function W_update = update(W, lambda, alpha, M, N, P, H_B, H_E, s, N_B, N_E)
    grad_W = gradient(W, lambda, M, H_B, H_E, s, N_B, N_E);

    W = W - alpha * grad_W;
    W_update = projection(W, N, P);
end

function grad_W = gradient(W, lambda, M, H_B, H_E, s, N_B, N_E)
    grad_W = -gradFunction(M, H_B, s, N_B, W) + lambda * gradFunction(M, H_E, s, N_E, W);
end

function W_projection = projection(W, N, P)
    w_vectorization = vectorizationCal(N, W);
    
    condition = norm(w_vectorization) / sqrt(P);

    if condition <= 1
        w_projection = w_vectorization / sqrt(P);
    else
        w_projection = w_vectorization / norm(w_vectorization);
    end

    W_projection = inverse_vectorizationCal(N, w_projection);
end

function vector = vectorizationCal(N, matrix)
    vector = reshape(matrix, N * N, 1);
end

function matrix = inverse_vectorizationCal(N, vector)
    matrix = reshape(vector, N, N);
end