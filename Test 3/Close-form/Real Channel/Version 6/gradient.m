function [grad_W] = gradient(W)
    global M; % M-ary modulation
    global H_B; % Channel matrix from Alice to Bob
    global H_E; % Channel matrix from Alice to Eve
    global s; % Modulated symbol vector
    global N_0; % Noise power
    global lambda; % Trade off parameter for objective function

    grad_W = -gradFunction(M,H_B,s,N_0,W) + lambda * gradFunction(M,H_E,s,N_0,W);
end