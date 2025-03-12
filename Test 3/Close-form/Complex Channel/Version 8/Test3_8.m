clear;
clc;

%==========================================================================
% Parameters
N_Comp = 2; % Number of the transmitted antennas
K_Comp = 2; % Number of the received antennas => Not depend on for this case! 

N = 2*N_Comp;
K = 2*K_Comp;
M = 4; % M-ary modulation
alpha = 0.0000025; % Step size
lambda = 0.5;
N_B = 0.000005; % Noise power
N_E = 0.000015; % Noise power
channel_var = 1;

num_iter = 1000; 
max_iteration = 50;
stop_criteria = 0.00001;

%==========================================================================
% Initial value 

W_0_Comp = complex(rand(N_Comp,N_Comp),rand(N_Comp,N_Comp)); % Generate randomly
W_0 = [real(W_0_Comp), -imag(W_0_Comp); imag(W_0_Comp), real(W_0_Comp)];

% s_Comp = [1 1 -1 -1; -1 1 -1 1]; % => BPSK modulation constellations {-1, 1}

s_Comp = [1+1i, -1-1i, -1+1i, -1-1i;...
          1-1i,  1-1i,  1-1i, -1+1i]; % => QPSK modulation constellations {+-1 +- 1i}

s = [real(s_Comp(:,1)) real(s_Comp(:,2)) real(s_Comp(:,3)) real(s_Comp(:,4));...
     imag(s_Comp(:,1)) imag(s_Comp(:,2)) imag(s_Comp(:,3)) imag(s_Comp(:,4))];

% Range of SNR values
SNR_dB = (-10):10:80; % dB scale
SNR = 10.^(SNR_dB/10); % Linear scale
P = SNR*N_B;

for iter = 1:num_iter 
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

    flag = 0;

    for j = 1:length(P)

        [W_update] = update(W,lambda,alpha,M,N,P(j),H_B,H_E,s,N_B,N_E); 
        
        objective_function_Bob(1,iter) = objectFunction(M,H_B,s,N_B,W_update); % Error Probability of Bob
        
        objective_function_Eve(1,iter) = objectFunction(M,H_E,s,N_E,W_update); % Error Probability of Eve
        
        objective_function(1,iter) = objective_function_Bob(1,iter) - lambda * objective_function_Eve(1,iter); % Objective Value
    
        objective_function_current = objective_function(1,iter);
        
        W = W_update;
    
        for i = 2:max_iteration 
    
            [W_update] = update(W,lambda,alpha,M,N,P(j),H_B,H_E,s,N_B,N_E); 
            
            objective_function_Bob(i,iter) = objectFunction(M,H_B,s,N_B,W_update); % Error Probability of Bob
        
            objective_function_Eve(i,iter) = objectFunction(M,H_E,s,N_E,W_update); % Error Probability of Eve
        
            objective_function(i,iter) = objective_function_Bob(i,iter) - lambda * objective_function_Eve(i,iter); % Objective Value
    
            if abs(objective_function(i,iter) - objective_function_current) <= stop_criteria
                if flag == 0 
                    objective_function_Bob_convergent(j,iter) = objective_function_Bob(i,iter);
                    objective_function_Eve_convergent(j,iter) = objective_function_Eve(i,iter);
                    flag = 1;
                end 
            end
            
            objective_function_current = objective_function(i,iter);
        
            W = W_update;
        end
    end
end

% Find mean of probability values
for j = 1:length(P)

    objective_function_Bob_convergent_average(j) = sum(objective_function_Bob_convergent(j,:))/num_iter;
    objective_function_Eve_convergent_average(j) = sum(objective_function_Eve_convergent(j,:))/num_iter;
end

%==========================================================================
%Figure of the objective value versus iteration 
figure(1)
semilogy(SNR_dB,objective_function_Bob_convergent_average,'r');
grid on
hold on
semilogy(SNR_dB,objective_function_Eve_convergent_average,'g');
hold off
xlabel('SNR');
ylabel('Error probability');
legend('Bob','Eve');