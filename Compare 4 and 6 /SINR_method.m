function [w_optimal] = SINR_method(H_B, H_E) 
    % =================================================================
    % Comparison: A. Mukherjee and A. L. Swindlehurst, "Robust Beamforming for Security in MIMO Wiretap Channels With Imperfect CSI," in IEEE Transactions on Signal Processing, vol. 59, no. 1, pp. 351-361, Jan. 2011.
          
    % For case: K_E >= N
    scaleFactor = 0.067;
            
    [eigenVectors,eigenValues] = eig(((H_B')*H_B),((H_E')*H_E));
            
    [max_value, max_index] = max(max(eigenValues));
            
    w_optimal = eigenVectors(:,max_index);
    
    if (norm(w_optimal) > 1) 
        w_optimal = w_optimal*scaleFactor;
    end
end