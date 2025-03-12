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

lambda = [1 0.5 1.5]; % Trade off parameter for objective function. It's noted that the gradient depends on lambda.

%==========================================================================
% Initial value 
max_iteration = 600;
stop_criteria = 0.0001;

objective_function_Bob=zeros(length(lambda),length(max_iteration));
objective_function_Eve=zeros(length(lambda),length(max_iteration));
objective_function=zeros(length(lambda),length(max_iteration));

for iter = 1:length(lambda)   
    % Optimization process
    W = W_0;

    % flag = 0;

    for i = 1:max_iteration 
        [W_update(:,:,iter)] = update(W,lambda(iter));
    
        objective_function_Bob(iter,i) = objectFunction(M,H_B,s,N_0,W_update(:,:,iter)); % Error Probability of Bob

        objective_function_Eve(iter,i) = objectFunction(M,H_E,s,N_0,W_update(:,:,iter)); % Error Probability of Eve

        % objective_function(iter,i) = objective_function_Bob(iter,i) - lambda(iter) * objective_function_Eve(iter,i); % Objective Value

        % if (norm(W_update(:,:,iter) - W,'fro')) <= stop_criteria
        %     if flag == 0 
        %         iter
        %         i
        %         W_update(:,:,iter)
        %         objective_function(iter,i)
        %         flag = 1;
        %     end 
        % end

        W = W_update(:,:,iter);
    end
end


%==========================================================================
%Figure of the objective value versus iteration 
figure(3)
iterations = 1:max_iteration;
plot(iterations,objective_function_Bob(1,:),'--b',iterations,objective_function_Eve(1,:),'-b');
grid on
hold on
plot(iterations,objective_function_Bob(2,:),'--r',iterations,objective_function_Eve(2,:),'-r');
plot(iterations,objective_function_Bob(3,:),'--g',iterations,objective_function_Eve(3,:),'-g');
hold off
xlabel('Iterations');
ylabel('Error probability');
title('Projected Gradient Descent');
legend("P_e^{Bob}, \lambda_1 = " + lambda(1),"P_e^{Eve}, \lambda_1 = " + lambda(1),...
       "P_e^{Bob}, \lambda_2 = " + lambda(2),"P_e^{Eve}, \lambda_2 = " + lambda(2),...
       "P_e^{Bob}, \lambda_3 = " + lambda(3),"P_e^{Eve}, \lambda_3 = " + lambda(3),'Location','northwest');