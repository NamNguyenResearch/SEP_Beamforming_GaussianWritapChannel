function [grad_W] = gradient(W,lambda,M,H_B,H_E,s,N_B,N_E)
    grad_W = -gradFunction(M,H_B,s,N_B,W) + lambda * gradFunction(M,H_E,s,N_E,W);
end