clear;
clc;

%==========================================================================
% Parameters
N_Comp = 2; % Number of the transmitted antennas
K_Comp = 2; % Number of the received antennas => Not depend on for this case! 

N = 2*N_Comp;
K = 2*K_Comp;
M = 4; % M-ary modulation
alpha = 0.005; % Step size
P = 1.5; % Given power
N_B = 0.005; % Noise power
N_E = 0.005; % Noise power
channel_var = 1;

num_iter = 5000; 
max_iteration = 25;

lambda = [0.5, 0.8, 1.2];

%==========================================================================
% Initial value 

W_0_Comp = complex(rand(N_Comp,N_Comp),rand(N_Comp,N_Comp)); % Generate randomly
W_0 = [real(W_0_Comp), -imag(W_0_Comp); imag(W_0_Comp), real(W_0_Comp)];

% s_Comp = [1 1 -1 -1; -1 1 -1 1]; % => BPSK modulation constellations {-1, 1}

s_Comp = [1+1i, -1-1i, -1+1i, -1-1i;...
          1-1i,  1-1i,  1-1i, -1+1i]; % => QPSK modulation constellations {+-1 +- 1i}

s = [real(s_Comp(:,1)) real(s_Comp(:,2)) real(s_Comp(:,3)) real(s_Comp(:,4));...
     imag(s_Comp(:,1)) imag(s_Comp(:,2)) imag(s_Comp(:,3)) imag(s_Comp(:,4))];

for iter = 1:num_iter
    % Generate the channel matrices and the modulated symbol vectors

    % Circularly symmetric complex Gaussian (Rayleigh fading channel realizations) channel m atrices
    H_B_Comp = sqrt(channel_var/2)*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp)); % Channel matrix from Alice to Bob
    H_E_Comp = sqrt(channel_var/2)*(randn(K_Comp,N_Comp)+1i*randn(K_Comp,N_Comp)); % Channel matrix from Alice to Eve            
            
    % % Scattering channel matrices
    % H_B_Comp = create_H(K_Comp,N_Comp,0.1,5,20); 
    % H_E_Comp = create_H(K_Comp,N_Comp,0.1,5,20);
            
    % Real channel matrices
    H_B = [real(H_B_Comp), -imag(H_B_Comp); imag(H_B_Comp), real(H_B_Comp)];
    H_E = [real(H_E_Comp), -imag(H_E_Comp); imag(H_E_Comp), real(H_E_Comp)];
            
    % Optimization process
    W = W_0;

    for i = 1:max_iteration

        for j = 1:length(lambda)  

            [W_update] = update(W,lambda(j),alpha,M,N,P,H_B,H_E,s,N_B,N_E); 
            
            objective_function_Bob(j,i,iter) = objectFunction(M,H_B,s,N_B,W_update); % Error Probability of Bob
        
            objective_function_Eve(j,i,iter) = objectFunction(M,H_E,s,N_E,W_update); % Error Probability of Eve
        
            W = W_update;
        end 
    end
end

% Find mean of probability values
for i = 1:max_iteration

    for j = 1:length(lambda)

        objective_function_Bob_average(j,i) = sum(objective_function_Bob(j,i,:))/num_iter;
        objective_function_Eve_average(j,i) = sum(objective_function_Eve(j,i,:))/num_iter;
    end
end

%==========================================================================
%Figure of the objective value versus iteration 
figure(1)
iterations = 1:max_iteration;
plot(iterations,objective_function_Bob_average(1,:),'-b',iterations,objective_function_Eve_average(1,:),'--b','LineWidth',1);
grid on
hold on
plot(iterations,objective_function_Bob_average(2,:),'-r',iterations,objective_function_Eve_average(2,:),'--r','LineWidth',1);
plot(iterations,objective_function_Bob_average(3,:),'-g',iterations,objective_function_Eve_average(3,:),'--g','LineWidth',1);
hold off
xlabel('Iterations');
ylabel('Bound of symbol error probability');
legend("P_e^{Bob}, \gamma_1 = " + lambda(1),"P_e^{Eve}, \gamma_1 = " + lambda(1),...
       "P_e^{Bob}, \gamma_2 = " + lambda(2),"P_e^{Eve}, \gamma_2 = " + lambda(2),...
       "P_e^{Bob}, \gamma_3 = " + lambda(3),"P_e^{Eve}, \gamma_3 = " + lambda(3),'Location','NorthWest');