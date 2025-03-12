function [objective_function_Bob_average, objective_function_Eve_average] = optimization_case(N_Comp,K_Comp,L_Comp)
    %==========================================================================
    % Parameters    
    N = 2*N_Comp;
    K = 2*K_Comp;
    L = 2*L_Comp;
   
    M = 4; % M-ary detection
    lambda = 1; % Tradeoff Parameter
    P = 1.5; % Given power
    N_B = 0.005; % Noise power
    N_E = 0.005; % Noise power
    channel_var = 1;
    
    num_iter = 5000; 
    max_iteration = 45;
    
    %==========================================================================
    % Initial value 
    s_Comp = [1+1i, -1-1i, -1+1i, -1-1i;...
              1-1i,  1-1i,  1-1i, -1+1i]; % => QPSK modulation constellations {+- 1 +- 1i}

    s = [real(s_Comp(:,1)) real(s_Comp(:,2)) real(s_Comp(:,3)) real(s_Comp(:,4));...
        imag(s_Comp(:,1)) imag(s_Comp(:,2)) imag(s_Comp(:,3)) imag(s_Comp(:,4))];
    for iter = 1:num_iter 
    
        W_0_Comp = complex(rand(N_Comp,L_Comp),rand(N_Comp,L_Comp)); % Generate randomly
    
        W_0 = [real(W_0_Comp), -imag(W_0_Comp); imag(W_0_Comp), real(W_0_Comp)];
    
        % Generate the channel matrices and the modulated symbol vectors
    
        % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel m atrices
        H_B_Comp = sqrt(channel_var/2)*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp));
        H_E_Comp = sqrt(channel_var/2)*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp));
        
        % % Scattering channel matrices
        % H_B_Comp = create_H(K_Comp,N_Comp,0.1,5,20); 
        % H_E_Comp = create_H(K_Comp,N_Comp,0.1,5,20);
        
        % Real channel matrices
        H_B = [real(H_B_Comp), -imag(H_B_Comp); imag(H_B_Comp), real(H_B_Comp)];
        H_E = [real(H_E_Comp), -imag(H_E_Comp); imag(H_E_Comp), real(H_E_Comp)];
        
        % Optimization process
        W = W_0;
        
        for i = 1:max_iteration 
    
            [W_update] = update(W,lambda,M,N,L,P,H_B,H_E,s,N_B,N_E); 
            
            objective_function_Bob(i,iter) = objectFunction(M,H_B,s,N_B,W_update); % Error Probability of Bob
        
            objective_function_Eve(i,iter) = objectFunction(M,H_E,s,N_E,W_update); % Error Probability of Eve
        
            W = W_update;
        end
    end
    
    % Find mean of probability values
    for i=1:max_iteration
    
        objective_function_Bob_average(i) = sum(objective_function_Bob(i,:))/num_iter;
        objective_function_Eve_average(i) = sum(objective_function_Eve(i,:))/num_iter;
    end
end
