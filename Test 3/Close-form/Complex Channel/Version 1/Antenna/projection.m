function [W_projection] = projection(W,N,L,P)

    w_vectorization = vectorizationCal(N,L,W);
    
    condition = norm(w_vectorization/(sqrt(P)));

    if (condition <= 1)
        w_projection = w_vectorization/(sqrt(P));
    else
        w_projection = w_vectorization/(norm(w_vectorization));
    
    end 

    W_projection = inverse_vectorizationCal(N,L,w_projection);
end