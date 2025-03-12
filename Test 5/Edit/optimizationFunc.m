function [W_optimal] = optimizationFunc(H_B,H_E,s_0,s_1,P,N,constraint)
     cvx_begin quiet
        variable W_optimal(N,N) complex
        
        minimize ( norm(H_E*W_optimal*(s_0-s_1)) )
        subject to 
            real(H_B*W_optimal*(s_0-s_1)) >= constraint
            (norm(W_optimal,'fro')) <= sqrt(P)
    cvx_end
end