clear;
clc;

%==========================================================================
% Parameters
global N;       % Number of the transmitted antennas
global K;       % Number of the received antennas => Not depend on for this case!
global M;       % M-ary modulation

global alpha;   % Step size
global P;       % Given power
global N_0;     % Noise power

global H_B;     % Channel matrix from Alice to Bob
global H_E;     % Channel matrix from Alice to Eve
global s;       % Modulated symbol vector

N = 2;
K = 2;
M = 4;

alpha = 0.005;
P = 1;
N_0 = 0.0001;

% Generate the channel matrices and the modulated symbol vectors
H_B = [0.1  0.02; 0.013  0.12]; % => Bob's channel is better than Eve's one 
H_E = [0.01 0.03; 0.015 0.001];

s = [1 1 -1 -1; -1 1 -1 1]; % => BPSK modulation constellations {-1, 1}

%==========================================================================
% Initial value 
% W_0 = randn(N);

W_0 = [0.5359  -0.6875;  0.6929, 0.3176];

lambda = 0.1:0.1:1.5; % Trade off parameter for objective function. It's noted that the gradient depends on lambda.

%==========================================================================
% Initial value 
max_iteration = 450;
stop_criteria = 0.0001;

for j = 1:length(lambda)   
    % Optimization process
    W = W_0;
    flag = 0;
    
    [W_update(:,:,j)] = update(W,lambda(j));
    
    objective_function_Bob(1,j) = objectFunction(M,H_B,s,N_0,W_update(:,:,j)); % Error Probability of Bob
    
    objective_function_Eve(1,j) = objectFunction(M,H_E,s,N_0,W_update(:,:,j)); % Error Probability of Eve
    
    objective_function(1,j) = objective_function_Bob(1,j) - lambda(j) * objective_function_Eve(1,j); % Objective Value
    
    objective_function_current = objective_function(1,j);
    
    W = W_update(:,:,j);

    for i = 2:max_iteration 
        [W_update(:,:,j)] = update(W,lambda(j));
    
        objective_function_Bob(i,j) = objectFunction(M,H_B,s,N_0,W_update(:,:,j)); % Error Probability of Bob
        
        objective_function_Eve(i,j) = objectFunction(M,H_E,s,N_0,W_update(:,:,j)); % Error Probability of Eve
        
        objective_function(i,j) = objective_function_Bob(i,j) - lambda(j) * objective_function_Eve(i,j); % Objective Value
        
        if flag == 0 
            if abs(objective_function(i,j) - objective_function_current) <= stop_criteria
                objective_function_convergent(j) = objective_function(i,j);
                objective_function_Bob_convergent(j) = objective_function_Bob(i,j);
                objective_function_Eve_convergent(j) = objective_function_Eve(i,j);
                flag = 1;
            end
        end

        objective_function_current = objective_function(i,j);

        W = W_update(:,:,j);
    end
end


%==========================================================================
figure(1)
plot(lambda,objective_function_Bob_convergent,lambda,objective_function_Eve_convergent);
grid on
xlabel('\lambda');
ylabel('Error probability');
legend('Bob','Eve');