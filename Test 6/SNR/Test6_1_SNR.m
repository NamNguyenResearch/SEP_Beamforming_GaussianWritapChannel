clear;
clc;

% Parameters
N = 2;
K = 2;
N_B = 0.01;
N_E = 0.01;
channel_var = 1;
a = 1;
P = 1;

D = 0.0001;
constraint = (sqrt(N_B/2)*qfuncinv(D)/abs(a)).^2;
num_iter = 100; % Number of Monte-Carlo runs

EbN0dB_vector = [-15 -12 -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10];
Eb = 1;

for iter = 1:num_iter
    % Generate the channel matrices and the modulated symbol vectors
    H_B = sqrt(channel_var/2) * (randn(K, N) + 1i * randn(K, N));
    H_E = sqrt(channel_var/2) * (randn(K, N) + 1i * randn(K, N));

    w_optimal = optimizationFunc(H_B, H_E, P, N, constraint);

    % Beamforming vector normalization
    w_optimal = w_optimal/norm(w_optimal);

    for snr_i = 1:length(EbN0dB_vector)
        EbN0dB = EbN0dB_vector(snr_i);
        EbN0 = 10^(EbN0dB / 10);
        N0 = Eb / EbN0;

        norm_term_Bob = norm(H_B * w_optimal * abs(a));
        norm_term_Eve = norm(H_E * w_optimal * abs(a));

        ErrorProbability_Bob(snr_i,iter) = qfunc(norm_term_Bob / sqrt(N0 / 2));
        ErrorProbability_Eve(snr_i,iter) = qfunc(norm_term_Eve / sqrt(N0 / 2));
    end
end

% Calculate averages
ErrorProbability_Bob_average = mean(ErrorProbability_Bob, 2);
ErrorProbability_Eve_average = mean(ErrorProbability_Eve, 2);

% Plotting
figure(1)
semilogy(EbN0dB_vector, ErrorProbability_Bob_average, 'b-o', EbN0dB_vector, ErrorProbability_Eve_average, 'r-*', 'LineWidth', 1.5);
grid on
xlabel('SNR (dB)');
ylabel('Symbol Error Probability');
legend('Bob','Eve','Location','SouthWest');