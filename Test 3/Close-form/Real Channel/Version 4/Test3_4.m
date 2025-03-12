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
number_test = 15;

for iter = 1: number_test 
    % Initial value 
    % iter
    W_0 = randn(N);

    max_iteration = 1000;
%     stop_criteria = 0.0001;

    % Optimization process
    W = W_0;
    objective_function=zeros(1,length(max_iteration));

    
%     flag = 0;

    for i = 1:max_iteration 
        [W_update] = update(W);

        objective_function(i) = objectFunction(M,H_B,s,N_0,W_update) ...
                              - lambda * objectFunction(M,H_E,s,N_0,W_update);

%         if (norm(W_update - W,'fro')) <= stop_criteria
%             if flag == 0 
%                 iter
%                 i
%                 W_update
%                 objective_function(i);
%                 flag = 1;
%             end 
%         end

        W = W_update;
    end


    %======================================================================
    %Figure of the objective value versus iteration 
    iterations = 1:max_iteration;
    plot(iterations,objective_function);
    grid on
    hold on;
    xlabel('Iterations');
    ylabel('Objective Value');
end
hold off