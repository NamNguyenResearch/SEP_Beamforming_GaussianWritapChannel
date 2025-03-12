function [W_projection] = projection(W)
    global N; % Number of the transmitted antennas
    global P; % Given power
    
    w_vectorization = vectorizationCal(N,W);
    
    condition = norm(w_vectorization/(sqrt(P)));

    if (condition <= 1)
        w_projection = w_vectorization/(sqrt(P));
    else
        w_projection = w_vectorization/(norm(w_vectorization));
    
    end 

    W_projection = inverse_vectorizationCal(N,w_projection);
end