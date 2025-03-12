clear;
clc;

%==========================================================================
% Parameters
global N;       % Number of the transmitted antennas
global K;       % Number of the received antennas => Not depend on for this case!
global M;       % M-ary modulation

global alpha;   % Step size
global lambda;  % Trade off parameter for objective function
global P;       % Given power
global N_0;     % Noise power

global H_B;     % Channel matrix from Alice to Bob
global H_E;     % Channel matrix from Alice to Eve
global s;       % Modulated symbol vector

N = 2;
K = 2;
M = 4;

alpha = 0.005;
lambda = 1;
P = 1;
N_0 = 0.0001;

% Generate the channel matrices and the modulated symbol vectors
H_B = [0.1  0.02; 0.013  0.12]; % => Bob's channel is better than Eve's one 
H_E = [0.01 0.03; 0.015 0.001];

s = [1 1 -1 -1; -1 1 -1 1]; % => BPSK modulation constellations {-1, 1}

%==========================================================================
number_test = 40; 
objective_function_Eve_constraint = 100000000;
max_iteration = 600;
stop_criteria = 0.00001;
j_min = 0;

j_convergent = zeros(1,length(number_test));

for j = 1:number_test 
    % Initial value 
    W_0 = randn(N);

    % Optimization process
    W = W_0;
    flag = 0;
    
    [W_update] = update(W);
    
    objective_function_Bob(1,j) = objectFunction(M,H_B,s,N_0,W_update); % Error Probability of Bob
    
    objective_function_Eve(1,j) = objectFunction(M,H_E,s,N_0,W_update); % Error Probability of Eve
    
    objective_function(1,j) = objective_function_Bob(1,j) - lambda * objective_function_Eve(1,j); % Objective Value
    
    objective_function_current = objective_function(1,j);
    
    W = W_update;

    for i = 2:max_iteration 
        [W_update] = update(W);
        
        objective_function_Bob(i,j) = objectFunction(M,H_B,s,N_0,W_update); % Error Probability of Bob

        objective_function_Eve(i,j) = objectFunction(M,H_E,s,N_0,W_update); % Error Probability of Eve

        objective_function(i,j) = objective_function_Bob(i,j) - lambda * objective_function_Eve(i,j); % Objective Value
        
        if flag == 0 
            if abs(objective_function(i,j) - objective_function_current) <= stop_criteria
                objective_function_convergent(j) = objective_function(i,j);
                objective_function_Eve_convergent(j) = objective_function_Eve(i,j);
                j_convergent(j) = i;
                flag = 1;
            end
        end
        
        objective_function_current = objective_function(i,j);

        W = W_update;
    end
    
    if (flag == 1) 
        if (objective_function_Eve_convergent(j) < objective_function_Eve_constraint) && (objective_function_Eve_convergent(j) <= 1) 
            j_min = j;
            objective_function_Eve_constraint = objective_function_Eve_convergent(j);
        end
    end 
end

%==========================================================================
%Figure of the objective value versus iteration 
iterations = 1:max_iteration;

plot(iterations,objective_function_Bob(:,j_min),'--r',iterations,objective_function_Eve(:,j_min),'-r');
hold on

for j = 1:number_test
    if (j_convergent(j) ~= 0) 
        if (j ~= j_min) && (objective_function_Eve_convergent(j) <= 1) 
            plot(iterations,objective_function_Bob(:,j),'--k',iterations,objective_function_Eve(:,j),'-k');
        end
    end
end

hold off
grid on;
xlabel('Iterations');
ylabel('Error Probability');
legend('Optimal Bob','Optimal Eve','Bob','Eve');