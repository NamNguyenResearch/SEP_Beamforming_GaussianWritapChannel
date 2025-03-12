function [projection_W] = projection(W)
    global N; % Number of the transmitted antennas
    global P; % Given power
    
    cvx_begin quiet
        variable projection_W(N,N) 
        
        minimize (norm(projection_W - W,'fro'))
        subject to 
            % trace(projection_W*(projection_W')) <= P  
            % trace(projection_W) <= P % Find W and use the decomposition
            % square_pos(norm(projection_W,'fro')) <= P

            (norm(projection_W,'fro')) <= sqrt(P)
    cvx_end
end