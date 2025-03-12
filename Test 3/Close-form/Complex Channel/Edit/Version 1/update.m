function W_update = update(W, lambda, alpha, M, N, P, H_B, H_E, s, N_B, N_E)
    grad_W = gradient(W, lambda, M, H_B, H_E, s, N_B, N_E);

    W = W - alpha * grad_W;
    W_update = projection(W, N, P);
end