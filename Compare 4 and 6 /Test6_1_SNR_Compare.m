clear;
clc;

% Parameters
N = 2;
K = 2;
N_B = 0.01;
N_E = 0.01;
channel_var = 1;
a = 1;

D_Bob = 0.0001;
P = 1;
constraint = (sqrt(N_B/2) * qfuncinv(D_Bob) / abs(a))^2;
num_iter = 1; % Number of Monte-Carlo runs

D_Eve = 0.3;
SNR_dB = 10;
SNR = 10^(SNR_dB / 10);
P2 = SNR * N_B;

EbN0dB_vector = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10];
Eb = 1;

num_SNR = length(EbN0dB_vector);

% Preallocate arrays for error probabilities
ErrorProbability_Bob = zeros(num_SNR, num_iter);
ErrorProbability_Eve = zeros(num_SNR, num_iter);
ErrorProbability_Bob2 = zeros(num_SNR, num_iter);
ErrorProbability_Eve2 = zeros(num_SNR, num_iter);
ErrorProbability_Bob3 = zeros(num_SNR, num_iter);
ErrorProbability_Eve3 = zeros(num_SNR, num_iter);

% Pre-generate channel matrices (assuming fixed for the given example)
H_B = [0.0262 0.0049; -0.1598 -0.2414]; 
H_E = [0.0498 0.0194; -0.0446 -0.0758]; 

for iter = 1:num_iter
    % Optimization functions (assuming these functions are defined elsewhere)
    w_optimal = optimizationFunc(H_B, H_E, P, N, constraint);
    w_optimal2 = eigen_decomposition_method(H_B, H_E, N, N_B, N_E, P2, D_Eve, a);
    w_optimal3 = SINR_method(H_B, H_E);

    % Beamforming vector normalization
    w_optimal = w_optimal / norm(w_optimal);
    w_optimal2 = w_optimal2 / norm(w_optimal2);
    w_optimal3 = w_optimal3 / norm(w_optimal3);
    
    for snr_i = 1:num_SNR
        EbN0dB = EbN0dB_vector(snr_i);
        EbN0 = 10^(EbN0dB / 10);
        N0 = Eb / EbN0;

        % Calculate norm terms
        norm_term_Bob = norm(H_B * w_optimal * abs(a));
        norm_term_Eve = norm(H_E * w_optimal * abs(a));

        norm_term_Bob2 = norm(H_B * w_optimal2 * abs(a));
        norm_term_Eve2 = norm(H_E * w_optimal2 * abs(a));

        norm_term_Bob3 = norm(H_B * w_optimal3 * abs(a));
        norm_term_Eve3 = norm(H_E * w_optimal3 * abs(a));
        
        % Calculate error probabilities
        ErrorProbability_Bob(snr_i, iter) = qfunc(norm_term_Bob / sqrt(N0 / 2));
        ErrorProbability_Eve(snr_i, iter) = qfunc(norm_term_Eve / sqrt(N0 / 2));

        ErrorProbability_Bob2(snr_i, iter) = qfunc(norm_term_Bob2 / sqrt(N0 / 2));
        ErrorProbability_Eve2(snr_i, iter) = qfunc(norm_term_Eve2 / sqrt(N0 / 2));

        ErrorProbability_Bob3(snr_i, iter) = qfunc(norm_term_Bob3 / sqrt(N0 / 2));
        ErrorProbability_Eve3(snr_i, iter) = qfunc(norm_term_Eve3 / sqrt(N0 / 2));
    end
end

% Calculate averages
ErrorProbability_Bob_average = mean(ErrorProbability_Bob, 2);
ErrorProbability_Eve_average = mean(ErrorProbability_Eve, 2);
ErrorProbability_Bob_average2 = mean(ErrorProbability_Bob2, 2);
ErrorProbability_Eve_average2 = mean(ErrorProbability_Eve2, 2);
ErrorProbability_Bob_average3 = mean(ErrorProbability_Bob3, 2);
ErrorProbability_Eve_average3 = mean(ErrorProbability_Eve3, 2);

% Plotting
figure(1)
plot(EbN0dB_vector, ErrorProbability_Bob_average, 'b-o', 'LineWidth', 1.5);
hold on
plot(EbN0dB_vector, ErrorProbability_Eve_average, 'r-o', 'LineWidth', 1.5);
plot(EbN0dB_vector, ErrorProbability_Bob_average2, 'b--', 'LineWidth', 1.5);
plot(EbN0dB_vector, ErrorProbability_Eve_average2, 'r--', 'LineWidth', 1.5);
plot(EbN0dB_vector, ErrorProbability_Bob_average3, 'b-*', 'LineWidth', 1.5);
plot(EbN0dB_vector, ErrorProbability_Eve_average3, 'r-*', 'LineWidth', 1.5);
hold off
grid on
xlabel('SNR (dB)');
ylabel('Symbol Error Probability');
legend('Bob, Egien-BF', 'Eve, Egien-BF', 'Bob, SDR-BF', 'Eve, SDR-BF', 'Bob, SINR-BF', 'Eve, SINR-BF', 'Location', 'SouthWest');