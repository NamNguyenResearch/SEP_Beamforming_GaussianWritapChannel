clear;
clc;

%==========================================================================
% Parameters
N = 2; % Number of transmitted antennas
K = 2; % Number of received antennas
M = 3; % M-ary modulation

alpha = 0.005; % Step size
lambda = 1;

N_B = 0.005; % Noise power
N_E = 0.005; % Noise power
channel_var = 1;

SNR_dB = -20:10:20;

SNR = 10.^(SNR_dB/10);

P = SNR * N_B;

num_iter = 2000;
max_iteration = 100;

%==========================================================================
% Initial value
W_0 = complex(rand(N, N), rand(N, N)); % Generate randomly

s = [1, -1, -1;...
     -1, 1, -1]; % => QPSK modulation constellations {+-1 +-1i}

% Preallocate arrays
objective_function_Bob = zeros(length(SNR_dB), num_iter);
objective_function_Eve = zeros(length(SNR_dB), num_iter);

for iter = 1:num_iter
    % Generate the channel matrices and the modulated symbol vectors
    H_B = sqrt(channel_var/2) * (randn(K, N) + 1i * randn(K, N));
    H_E = sqrt(channel_var/2) * (randn(K, N) + 1i * randn(K, N));

    for i = 1:length(SNR_dB)
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