function w_optimal = optimizationFunc(H_B, H_E, P, N, constraint)
    % Solve the semidefinite relaxation problem
    C_1 = ctranspose(H_E) * H_E;
    C_2 = ctranspose(H_B) * H_B;

    cvx_begin quiet
        variable W_optimal(N, N) hermitian

        minimize(trace(C_1 * W_optimal));
        subject to
            trace(C_2 * W_optimal) >= constraint;
            trace(W_optimal) <= P;
            W_optimal == semidefinite(N);
    cvx_end

    % Find the optimal beamforming vector
    w_optimal = randA(W_optimal, N, H_E, H_B, constraint, P);

    % Uncomment the desired method from below
    % w_optimal = randB(W_optimal, N, H_E, H_B, constraint, P);
    % w_optimal = randC(W_optimal, N, H_E, H_B, constraint, P);
    % w_optimal = rand_GaussianVector(W_optimal, N, H_E, H_B, constraint, P);
end