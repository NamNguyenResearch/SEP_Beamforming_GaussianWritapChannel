clear;
clc;

%==========================================================================
% Parameters
N_Comp = 2; % Number of transmitted antennas
K_Comp = 2; % Number of received antennas

N = 2 * N_Comp;
K = 2 * K_Comp;
M = 4; % M-ary modulation

alpha = 0.005; % Step size
lambda = 20; % Tradeoff Parameter

N_B = 0.005; % Noise power
N_E = 0.005; % Noise power
channel_var = 1;

SNR_dB = -30:10:10;
SNR = 10.^(SNR_dB/10);

P = SNR*N_B;

num_iter = 2000; 
max_iteration = 50;

%==========================================================================
% Initial value
W_0_Comp = complex(rand(N_Comp, N_Comp), rand(N_Comp, N_Comp)); % Generate randomly
W_0 = [real(W_0_Comp), -imag(W_0_Comp); imag(W_0_Comp), real(W_0_Comp)];

% s_Comp = [1+1i, -1-1i, -1+1i, -1-1i;...
%           1-1i,  1-1i,  1-1i, -1+1i]; % => QPSK modulation constellations {+-1 +- 1i}

s = [real(s_Comp(:, 1)) real(s_Comp(:, 2)) real(s_Comp(:, 3)) real(s_Comp(:, 4));...
     imag(s_Comp(:, 1)) imag(s_Comp(:, 2)) imag(s_Comp(:, 3)) imag(s_Comp(:, 4))];

% Preallocate arrays
objective_function_Bob = zeros(length(SNR_dB), num_iter);
objective_function_Eve = zeros(length(SNR_dB), num_iter);

for i = 1:length(SNR_dB)
    % Optimization process
    W = W_0;
    
    for iter = 1:num_iter
        % Generate the channel matrices and the modulated symbol vectors
        H_B_Comp = sqrt(channel_var/2) * (randn(K_Comp, N_Comp) + 1i * randn(K_Comp, N_Comp));
        H_E_Comp = sqrt(channel_var/2) * (randn(K_Comp, N_Comp) + 1i * randn(K_Comp, N_Comp));
        H_B = [real(H_B_Comp), -imag(H_B_Comp); imag(H_B_Comp), real(H_B_Comp)];
        H_E = [real(H_E_Comp), -imag(H_E_Comp); imag(H_E_Comp), real(H_E_Comp)];

        % Optimization process
        W = W_0;

        for j = 1:max_iteration
            W_update = update(W, lambda, alpha, M, N, P(i), H_B, H_E, s, N_B, N_E);

            objective_function_Bob(i, iter) = objectFunction(M, H_B, s, N_B, W_update);
            objective_function_Eve(i, iter) = objectFunction(M, H_E, s, N_E, W_update);

            W = W_update;
        end
    end
end

% Find mean of probability values
objective_function_Bob_average = mean(objective_function_Bob, 2);
objective_function_Eve_average = mean(objective_function_Eve, 2);

%==========================================================================
% Figure of the objective value versus iteration
figure(1);
semilogy(SNR_dB, objective_function_Bob_average, '-b','LineWidth',1.5);
hold on
semilogy(SNR_dB, objective_function_Eve_average, '-r','LineWidth',1.5);
hold off
grid on
xlabel('SNR');
ylabel('Error probability');
legend('Bob','Eve','Location','SouthEast');